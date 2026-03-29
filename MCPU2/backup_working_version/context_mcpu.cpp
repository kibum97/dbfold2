#include "context_mcpu.h"
#include "protein_structure.h"
#include "amino_acids.h"
#include "utils/geometry_utils.h"
#include "utils/mcpu_hbond_utils.h"
#include <iostream>
#include <chrono>

ContextMCPU::ContextMCPU(SystemMCPU& system, const Eigen::Matrix3Xf positions)
    : system_(system), positions_(positions), cell_list_(CellList(6.0f)) {
    // Resize vectors per chain
    n_chains_ = system_.n_residues_per_chain.size();
    backbone_plane_vectors_per_chain.resize(n_chains_);
    backbone_bisect_vectors_per_chain.resize(n_chains_);
    plane_vectors_per_chain.resize(n_chains_);
    bisect_vectors_per_chain.resize(n_chains_);
    backbone_phi_angles_per_chain.resize(n_chains_);
    backbone_psi_angles_per_chain.resize(n_chains_);
    backbone_phi_feats_per_chain.resize(n_chains_);
    backbone_psi_feats_per_chain.resize(n_chains_);
    backbone_plane_feats_per_chain.resize(n_chains_);
    backbone_bisect_feats_per_chain.resize(n_chains_);
    sidechain_torsion_feats_per_chain.resize(n_chains_);
    n_residues_ = 0;
    n_atoms_ = 0;
    for (int i = 0; i < n_chains_; ++i) {
        n_residues_ += system_.n_residues_per_chain[i];
        n_atoms_ += system_.n_atoms_per_chain[i];
    }
    // Initialize matrices
    pairwise_contact_matrix = Eigen::MatrixXf::Zero(n_atoms_, n_atoms_);
    ca_ca_distance2_matrix = Eigen::MatrixXf::Zero(n_residues_, n_residues_);
    hydrogen_bond_feature_matrix = Eigen::MatrixXi::Zero(n_residues_, n_residues_);

    // Constructor implementation
    reorder_positions();
    cell_list_.build_cells(positions_);
    // Initialize context per chain
    compute_pairwise_contacts();
    compute_backbone_vectors();
    compute_backbone_angles();
    compute_sidechain_angles();
    compute_aromatic_centers();
    compute_hydrogen_bonds();
}

ContextMCPU::~ContextMCPU() {
    // Destructor implementation
}

void ContextMCPU::initialize() {
    // Initialization code
}

void ContextMCPU::update() {
    // Update code
}

void ContextMCPU::reorder_positions() {
    // Reorder positions based on the system's atom indices
    Eigen::Matrix3Xf new_positions(3, system_.mcpu_atoms.size());
    for (size_t i = 0; i < system_.mcpu_atoms.size(); ++i) {
        int atom_id = system_.mcpu_atoms[i].atom->id;
        new_positions.col(i) = positions_.col(atom_id);
    }
    this->positions_ = new_positions;
}

void ContextMCPU::compute_pairwise_contacts() {
    for (int atom_id = 0; atom_id < positions_.cols(); ++atom_id) {
        mcpu_atom& atom = system_.mcpu_atoms[atom_id];
        Eigen::Vector3f atom_coords = positions_.col(atom_id);
        Eigen::Vector3i cell_index = cell_list_.get_cell_index(atom_coords);
        std::vector<int> neighbor_atom_ids = cell_list_.get_neighbor_atom_indices(atom_coords);
        bool is_backbone = (atom.smog_type > 78); // smog_type > 78 indicates a backbone atom
        for (int neighbor_atom_id : neighbor_atom_ids) {
            if (neighbor_atom_id != atom_id) {
                mcpu_atom& neighbor_atom = system_.mcpu_atoms[neighbor_atom_id];
                if (neighbor_atom.chain_id == atom.chain_id && abs(neighbor_atom.residue_id - atom.residue_id) < 4) {
                    continue; // Skip atoms within the same residue or close residues
                } else if (neighbor_atom.smog_type > 78 && is_backbone) {
                    continue; // Skip interactions between backbone atoms
                } else {
                    Eigen::Vector3f neighbor_coords = positions_.col(neighbor_atom_id);
                    float cutoff2 = system_.contact_cutoff2(atom.smog_type, neighbor_atom.smog_type);
                    float distance2 = (atom_coords - neighbor_coords).squaredNorm();
                    if (distance2 < cutoff2) {
                        if (atom_id > neighbor_atom_id) {
                            pairwise_contact_matrix(atom_id, neighbor_atom_id) = 1.0f;
                        } else {
                            pairwise_contact_matrix(neighbor_atom_id, atom_id) = 1.0f;
                        }
                    }
                    if (atom.atom->name == "CA" && neighbor_atom.atom->name == "CA") {
                        int residue_id = atom.residue_id;
                        int neighbor_residue_id = neighbor_atom.residue_id;
                        if (residue_id > neighbor_residue_id) {
                            std::swap(residue_id, neighbor_residue_id);
                        }
                        ca_ca_distance2_matrix(residue_id, neighbor_residue_id) = distance2;
                    }
                }
            }
        }
    }
}


void ContextMCPU::compute_backbone_vectors() {
    // Backbone atom indices are stored in system_.backbone_atoms
    // N will be 4 * i, CA will be 4 * i + 1, C will be 4 * i + 2, and O will be 4 * i + 3
    for (int chain_id = 0; chain_id < system_.n_residues_per_chain.size(); ++chain_id) {
        int n_residues = system_.n_residues_per_chain[chain_id];
        backbone_plane_vectors_per_chain[chain_id].resize(n_residues);
        backbone_bisect_vectors_per_chain[chain_id].resize(n_residues);
        plane_vectors_per_chain[chain_id].resize(3 * (n_residues - 1));
        bisect_vectors_per_chain[chain_id].resize(3 * (n_residues - 1));
        auto& backbone_atoms = system_.backbone_atoms[chain_id];
        for (int residue_id = 0; residue_id < n_residues; ++residue_id) {
            backbone_plane_vectors_per_chain[chain_id][residue_id] = compute_plane_vector(
                positions_.col(backbone_atoms[residue_id * 4]),
                positions_.col(backbone_atoms[residue_id * 4 + 1]),
                positions_.col(backbone_atoms[residue_id * 4 + 3])
            );
            backbone_bisect_vectors_per_chain[chain_id][residue_id] = compute_bisect_vector(
                positions_.col(backbone_atoms[residue_id * 4]),
                positions_.col(backbone_atoms[residue_id * 4 + 1]),
                positions_.col(backbone_atoms[residue_id * 4 + 3])
            );
        }
        for (int residue_id = 0; residue_id < n_residues - 1; ++residue_id) {
            plane_vectors_per_chain[chain_id][residue_id * 3] = compute_plane_vector(
                positions_.col(backbone_atoms[residue_id * 4]),
                positions_.col(backbone_atoms[residue_id * 4 + 1]),
                positions_.col(backbone_atoms[residue_id * 4 + 2])
            );
            bisect_vectors_per_chain[chain_id][residue_id * 3] = compute_bisect_vector(
                positions_.col(backbone_atoms[residue_id * 4]),
                positions_.col(backbone_atoms[residue_id * 4 + 1]),
                positions_.col(backbone_atoms[residue_id * 4 + 2])
            );
            plane_vectors_per_chain[chain_id][residue_id * 3  + 1] = compute_plane_vector(
                positions_.col(backbone_atoms[residue_id * 4 + 1]),
                positions_.col(backbone_atoms[residue_id * 4 + 2]),
                positions_.col(backbone_atoms[residue_id * 4 + 4])
            );
            bisect_vectors_per_chain[chain_id][residue_id * 3  + 1] = compute_bisect_vector(
                positions_.col(backbone_atoms[residue_id * 4 + 1]),
                positions_.col(backbone_atoms[residue_id * 4 + 2]),
                positions_.col(backbone_atoms[residue_id * 4 + 4])
            );
            plane_vectors_per_chain[chain_id][residue_id * 3  + 2] = compute_plane_vector(
                positions_.col(backbone_atoms[residue_id * 4 + 2]),
                positions_.col(backbone_atoms[residue_id * 4 + 4]),
                positions_.col(backbone_atoms[residue_id * 4 + 5])
            );
            bisect_vectors_per_chain[chain_id][residue_id * 3  + 2] = compute_bisect_vector(
                positions_.col(backbone_atoms[residue_id * 4 + 2]),
                positions_.col(backbone_atoms[residue_id * 4 + 4]),
                positions_.col(backbone_atoms[residue_id * 4 + 5])
            );
        }
    }
}

void ContextMCPU::compute_backbone_angles() {
    // Interate per chain
    for (int chain_id = 0; chain_id < system_.n_residues_per_chain.size(); ++chain_id) {
        // Resize the feature vectors
        backbone_phi_angles_per_chain[chain_id].resize(system_.n_residues_per_chain[chain_id] - 2);
        backbone_psi_angles_per_chain[chain_id].resize(system_.n_residues_per_chain[chain_id] - 2);
        backbone_phi_feats_per_chain[chain_id].resize(system_.n_residues_per_chain[chain_id] - 2);
        backbone_psi_feats_per_chain[chain_id].resize(system_.n_residues_per_chain[chain_id] - 2);
        backbone_plane_feats_per_chain[chain_id].resize(system_.n_residues_per_chain[chain_id] - 2);
        backbone_bisect_feats_per_chain[chain_id].resize(system_.n_residues_per_chain[chain_id] - 2);
        // Get backbone atom indices for the chain
        auto& backbone_atoms = system_.backbone_atoms[chain_id];
        // Iterate over the residues
        for (size_t i = 1; i < system_.n_residues_per_chain[chain_id] - 1; ++i) {
            // Get the backbone atoms
            Eigen::Vector3f n = positions_.col(backbone_atoms[i * 4]);
            Eigen::Vector3f ca = positions_.col(backbone_atoms[i * 4 + 1]);
            Eigen::Vector3f c = positions_.col(backbone_atoms[i * 4 + 2]);
            // Compute the angles
            // Phi angle and convert to index
            Eigen::Vector3f prev_c = positions_.col(backbone_atoms[(i - 1) * 4 + 2]);
            float phi_angle = compute_dihedral_angle(prev_c, n, ca, c);
            phi_angle = (phi_angle + M_PI) * (180.0f / M_PI);
            backbone_phi_angles_per_chain[chain_id][i - 1] = phi_angle;
            backbone_phi_feats_per_chain[chain_id][i - 1] = static_cast<int>(phi_angle/60.0f);
            // Psi angle and convert to index
            Eigen::Vector3f next_n = positions_.col(backbone_atoms[(i + 1) * 4]);
            float psi_angle = compute_dihedral_angle(n, ca, c, next_n);
            psi_angle = (psi_angle + M_PI) * (180.0f / M_PI);
            backbone_psi_angles_per_chain[chain_id][i - 1] = psi_angle;
            backbone_psi_feats_per_chain[chain_id][i - 1] = static_cast<int>(psi_angle/60.0f);
            // Plane angle and convert to index
            float plane_angle = compute_angle(
                backbone_plane_vectors_per_chain[chain_id][i - 1],
                backbone_plane_vectors_per_chain[chain_id][i + 1]
            );
            plane_angle = plane_angle * (180.0f / M_PI);
            backbone_plane_feats_per_chain[chain_id][i - 1] = static_cast<int>(plane_angle/30.0f);
            // Bisect angle and convert to index
            float bisect_angle = compute_angle(
                backbone_bisect_vectors_per_chain[chain_id][i - 1],
                backbone_bisect_vectors_per_chain[chain_id][i + 1]
            );
            bisect_angle = bisect_angle * (180.0f / M_PI);
            backbone_bisect_feats_per_chain[chain_id][i - 1] = static_cast<int>(bisect_angle/30.0f);
        }
    }
}

void ContextMCPU::compute_sidechain_angles() {
    // Compute side chain angles
    for (int chain_id = 0; chain_id < system_.n_residues_per_chain.size(); ++chain_id) {
        auto& backbone_atoms = system_.backbone_atoms[chain_id];
        // Resize the feature vectors
        int n_residues = system_.n_residues_per_chain[chain_id];
        sidechain_torsion_feats_per_chain[chain_id].resize(n_residues - 2);
        for (size_t i = 1; i < n_residues - 1; ++i) {
            // Get the backbone atoms
            int n_dihedrals = 0;
            for (n_dihedrals; n_dihedrals < system_.number_of_sidechain_dihedrals_per_chain[chain_id][i]; ++n_dihedrals) {
                // Compute the dihedral angle
                float chi_angle = compute_dihedral_angle(
                    positions_.col(backbone_atoms[i * 4] + n_dihedrals),
                    positions_.col(backbone_atoms[i * 4] + n_dihedrals + 1),
                    positions_.col(backbone_atoms[i * 4] + n_dihedrals + 2),
                    positions_.col(backbone_atoms[i * 4] + n_dihedrals + 3)
                );
                sidechain_torsion_feats_per_chain[chain_id][i - 1][n_dihedrals] = angle2index(chi_angle, 30.0f);
            }
            for (int total_dihedrals = 0; total_dihedrals < 4 - n_dihedrals; ++total_dihedrals) {
                sidechain_torsion_feats_per_chain[chain_id][i - 1][n_dihedrals + total_dihedrals] = 0;
            }
        }
    }
}

void ContextMCPU::compute_aromatic_centers() {
    // Compute the center of the aromatic ring
    for (const auto& aromatic_atoms : system_.aromatic_interactions) {
        const std::array<int, 3>& aromatic_atom_ids = aromatic_atoms.aromatic_atom_ids;
        aromatic_centers.push_back(compute_center({
            positions_.col(aromatic_atom_ids[0]),
            positions_.col(aromatic_atom_ids[1]),
            positions_.col(aromatic_atom_ids[2])
        }));   
    }
}

void ContextMCPU::compute_aromatic_angles() {
    for (auto& aromatic_center_1 : aromatic_centers) {
        for (auto& aromatic_center_2 : aromatic_centers) {
            if (&aromatic_center_1 != &aromatic_center_2) {
                // Compute the angle between the two aromatic centers
                float angle = compute_angle(aromatic_center_1, aromatic_center_2);
                // Convert to index
                if (angle > 89.9) {
                    angle = 89.9;
                }
                int aromatic_index = static_cast<int>(angle/10.);
                aromatic_interaction_feats.push_back(aromatic_index);
            }
        }
    }
}

void ContextMCPU::compute_hydrogen_bonds() {
    // Check if hydrogen bonds are formed
    for (int chain_id1 = 0; chain_id1 < system_.n_residues_per_chain.size(); ++chain_id1) {
        for (int chain_id2 = 0; chain_id2 < system_.n_residues_per_chain.size(); ++chain_id2) {
            for (int residue_id1 = 1; residue_id1 < system_.n_residues_per_chain[chain_id1] - 1; ++residue_id1) {
                for (int residue_id2 = 1; residue_id2 < system_.n_residues_per_chain[chain_id2] - 1; ++residue_id2) {
                    // Check if same chain
                    int residue_diff = get_residue_difference(chain_id1, residue_id1, chain_id2, residue_id2);
                    int donor_global_residue_id = system_.get_global_residue_id(chain_id1, residue_id1);
                    int acceptor_global_residue_id = system_.get_global_residue_id(chain_id2, residue_id2);
                    if (residue_diff < 4) {
                        hydrogen_bond_feature_matrix(donor_global_residue_id, acceptor_global_residue_id) = -1;
                        continue; // Skip if residues are too close
                    }
                    if (check_hydrogen_bond(chain_id1, residue_id1, acceptor_global_residue_id, chain_id2, residue_id2, acceptor_global_residue_id,
                                            backbone_phi_angles_per_chain, backbone_psi_angles_per_chain, ca_ca_distance2_matrix)) {
                        std::array<int, 3> bisect_indices = compute_angular_features(chain_id1, residue_id1, chain_id2, residue_id2, bisect_vectors_per_chain);
                        std::array<int, 3> plane_indices = compute_angular_features(chain_id1, residue_id1, chain_id2, residue_id2, plane_vectors_per_chain);
                        int secondary_structure = 0;
                        if (residue_diff != 4) {
                            secondary_structure = compute_secondary_structure(
                                positions_.col(system_.backbone_atoms[chain_id1][(residue_id1 + 1) * 3 + 1])
                                - positions_.col(system_.backbone_atoms[chain_id1][(residue_id1 - 1) * 3 + 1]),
                                positions_.col(system_.backbone_atoms[chain_id2][(residue_id2 + 1) * 3 + 1])
                                - positions_.col(system_.backbone_atoms[chain_id2][(residue_id2 - 1) * 3 + 1])
                            );
                        }
                        int hbond_index = compute_hbond_index(secondary_structure, bisect_indices, plane_indices);
                        hydrogen_bond_feature_matrix(donor_global_residue_id, acceptor_global_residue_id) = hbond_index;
                    } else {
                        hydrogen_bond_feature_matrix(donor_global_residue_id, acceptor_global_residue_id) = -1;
                    }
                }
            }
        }
    }
}

void ContextMCPU::compute_energy() {
    // Compute Mu potential energy
    auto mu_start = std::chrono::high_resolution_clock::now();
    Eigen::MatrixXf result = system_.mu_potential_coefficients.array() * pairwise_contact_matrix.array();
    energy.mu_energy = result.sum();
    auto mu_end = std::chrono::high_resolution_clock::now();
    // Compute Backbone Torsion energy
    float torsion_energy = 0.0;
    for (int chain_id = 0; chain_id < n_chains_; ++chain_id) {
        for (int res_id = 1; res_id < system_.n_residues_per_chain[chain_id] - 1; ++res_id) {
            int phi_index = backbone_phi_feats_per_chain[chain_id][res_id - 1];
            int psi_index = backbone_psi_feats_per_chain[chain_id][res_id - 1];
            int plane_index = backbone_plane_feats_per_chain[chain_id][res_id - 1];
            int bisect_index = backbone_bisect_feats_per_chain[chain_id][res_id - 1];  
            torsion_energy += system_.backbone_torsion_parameters_per_chain[chain_id][res_id - 1][plane_index][bisect_index][phi_index][psi_index];
        }
    }
    torsion_energy = torsion_energy/1000.; 
    energy.bb_torsion_energy = torsion_energy;
    auto bb_torsion_end = std::chrono::high_resolution_clock::now();
    // Compute Sidechain torsion energy
    float sc_torsion_energy = 0.0;
    for (int chain_id = 0; chain_id < n_chains_; ++chain_id) {
        for (int res_id = 1; res_id < system_.n_residues_per_chain[chain_id] - 1; ++res_id) {
            std::array<int, 4> chi_indices = sidechain_torsion_feats_per_chain[chain_id][res_id];
            std::cout << "DEBUG: Chain " << chain_id << ", Residue " << res_id
                      << " : "
                      << system_.residue_types_per_chain[chain_id][res_id - 1]
                      << " - "
                      << system_.residue_types_per_chain[chain_id][res_id]
                      << " - "
                      << system_.residue_types_per_chain[chain_id][res_id + 1]
                      << ": Chi indices = ["
                      << chi_indices[0] << ", "
                      << chi_indices[1] << ", "
                      << chi_indices[2] << ", "
                      << chi_indices[3] << "]"
                      << ", Energy = " << system_.sidechain_torsion_parameters_per_chain[chain_id][res_id][chi_indices[0]][chi_indices[1]][chi_indices[2]][chi_indices[3]]
                      << std::endl;
            sc_torsion_energy += system_.sidechain_torsion_parameters_per_chain[chain_id][res_id][chi_indices[0]][chi_indices[1]][chi_indices[2]][chi_indices[3]];
        }
    }
    sc_torsion_energy = sc_torsion_energy/1000.;
    energy.sc_torsion_energy = sc_torsion_energy;
    auto sc_torsion_end = std::chrono::high_resolution_clock::now();
    std::cout << "DEBUG: Sidechain torsion energy computed: " 
              << energy.sc_torsion_energy << std::endl;
    // Compute Hydrogen bond energy
    float hbond_energy = 0.0;
    // HBond - sequence specific
    for (int chain_id1 = 0; chain_id1 < system_.n_residues_per_chain.size(); ++chain_id1) {
        for (int chain_id2 = 0; chain_id2 < system_.n_residues_per_chain.size(); ++chain_id2) {
            for (int residue_id1 = 1; residue_id1 < system_.n_residues_per_chain[chain_id1] - 1; ++residue_id1) {
                for (int residue_id2 = 1; residue_id2 < system_.n_residues_per_chain[chain_id2] - 1; ++residue_id2) {
                    int i = system_.get_global_residue_id(chain_id1, residue_id1);
                    int j = system_.get_global_residue_id(chain_id2, residue_id2);
                    int hbond_index = hydrogen_bond_feature_matrix(i, j);
                    if (hbond_index >= 0) {
                        int secondary_structure_label = hbond_index / 531441; // 9^6 = 531441
                        if (system_.hydrogen_bond_parameters.find(hbond_index) == system_.hydrogen_bond_parameters.end()) {
                            continue; // Skip if the index is not found
                        }
                        hbond_energy += system_.sequence_based_hydorgen_bond_parameters[secondary_structure_label][system_.residue_types_per_chain[chain_id1][residue_id1]][system_.residue_types_per_chain[chain_id2][residue_id2]]
                                        * system_.hydrogen_bond_parameters.at(hbond_index);
                    }
                }
            }
        }
    }
    hbond_energy = hbond_energy/1000. * 2.;
    energy.hbond_energy = hbond_energy;
    auto hbond_end = std::chrono::high_resolution_clock::now();
    std::cout << "DEBUG: Hydrogen bond energy computed: " 
              << energy.hbond_energy << std::endl;
    // Compute Aromatic energy
    float aromatic_energy = 0.0;
    for (auto aromatic_feature : aromatic_interaction_feats) {
        aromatic_energy += system_.aromatic_parameters.at(aromatic_feature);
    }
    aromatic_energy = aromatic_energy/1000.;
    energy.aromatic_energy = aromatic_energy;
    auto aromatic_end = std::chrono::high_resolution_clock::now();
    std::cout << "DEBUG: Aromatic energy computed: " 
              << energy.aromatic_energy << std::endl;
    // Compute total energy
    energy.total_energy = energy.mu_energy + energy.bb_torsion_energy + energy.sc_torsion_energy + hbond_energy + aromatic_energy;

    // Calculate elapsed time for each part
    std::chrono::duration<double> mu_elapsed = mu_end - mu_start;
    std::chrono::duration<double> bb_torsion_elapsed = bb_torsion_end - mu_end;
    std::chrono::duration<double> sc_torsion_elapsed = sc_torsion_end - bb_torsion_end;
    std::chrono::duration<double> hbond_elapsed = hbond_end - sc_torsion_end;
    std::chrono::duration<double> aromatic_elapsed = aromatic_end - hbond_end;

    // Print elapsed time for each part
    std::cout << "Elapsed time for Mu potential: " << mu_elapsed.count() << " seconds" << std::endl;
    std::cout << "Elapsed time for Backbone Torsion: " << bb_torsion_elapsed.count() << " seconds" << std::endl;
    std::cout << "Elapsed time for Sidechain Torsion: " << sc_torsion_elapsed.count() << " seconds" << std::endl;
    std::cout << "Elapsed time for HBond: " << hbond_elapsed.count() << " seconds" << std::endl;
    std::cout << "Elapsed time for Aromatic: " << aromatic_elapsed.count() << " seconds" << std::endl;
}