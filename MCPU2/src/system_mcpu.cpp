#include "read_config.h"
#include "system_mcpu.h"
#include "protein_structure.h"
#include "amino_acids.h"
#include <unordered_map>
#include <iostream>


SystemMCPU::SystemMCPU(ForcefieldMCPU& forcefield, ProteinStructure& protein_structure)
    : forcefield_(forcefield), protein_structure_(protein_structure)
{
    // Constructor implementation
    n_chains_ = protein_structure_.num_chains;
    n_residues_ = protein_structure_.num_residues;
    n_atoms_ = protein_structure_.num_atoms;
    n_residues_per_chain = protein_structure_.num_residues_per_chain;
    n_atoms_per_chain = protein_structure_.num_atoms_per_chain;
    contact_cutoff2 = forcefield_.contact_cutoff2;

    // Initialize the real atoms
    pdb_atoms_to_mcpu_atoms();

    for (const auto& chain: protein_structure_.chains) {
        for (const auto& residue : chain->residues) {
            int chain_id = chain->id;
            int residue_id = residue->id;
            int global_residue_index = get_global_residue_index(chain_id, residue_id) + 1;
        }
    }
    init_backbone_atoms();
    init_sidechain_dihedrals();
    init_aromatic_interactions();
    init_hydrogen_bonds();
    init_system_forces();
}
SystemMCPU::~SystemMCPU() {
    // Destructor implementation
}

void SystemMCPU::pdb_atoms_to_mcpu_atoms() {
    // Convert PDB atoms to real atoms
    mcpu_atoms.clear();
    mcpu_atoms.resize(n_atoms_);
    chain_indices.resize(n_chains_);
    residue_indices.resize(n_residues_);
    std::string sidechain_key, backbone_key;
    std::unordered_map<std::string, std::pair<int, float>> smogTypeMap_ = forcefield_.smogTypeMap;
    for (const auto& chain : protein_structure_.chains) {
        int chain_id = chain->id;
        for (const auto& residue : chain->residues) {
            int residue_id = residue->id;
            std::unordered_map<std::string, int> sidechain_atom_offset = sidechain_atom_offset_map.at(residue->name);
            int global_residue_index = get_global_residue_index(chain_id, residue_id);
            int global_residue_id = get_global_residue_id(chain_id, residue_id);
            // Reorder the atoms in the residues to MCPU format
            // N -> CA -> sidechain atoms -> C -> O -> (OXT)
            // First count the backbone and sidechain atoms
            // Store atom to vector
            for (const auto& atom : residue->atoms) {
                // Get Smog Type and Radius
                // TODO: Instead of using XXX for backbone create a separate map for backbone atoms
                if (residue->name == "GLY" && atom->name == "CA") {
                    // Special case for Glycine, which has no sidechain
                    backbone_key = "GLY_CA";
                    sidechain_key = "GLY_NONE"; // This should not exist
                } else {
                    backbone_key = std::string("XXX_") + atom->name;
                    sidechain_key = residue->name + std::string("_") + atom->name;
                }
                // Current n_atoms_ and n_residues_ are used to assign the atom ID and the Residue ID
                // This is to give unique IDs to atoms and residues in the multi chain system
                if (smogTypeMap_.find(sidechain_key) != smogTypeMap_.end()) {
                    int smogType = smogTypeMap_[sidechain_key].first;
                    float smogRadius = smogTypeMap_[sidechain_key].second;
                    // Store at correct index in the mcpu_atoms vector
                    int sidechain_offset = sidechain_atom_offset[atom->name];
                    mcpu_atom current_atom = mcpu_atom(chain_id, global_residue_id, global_residue_index + sidechain_offset + 2, atom.get(), residue.get(), smogType, smogRadius);
                    mcpu_atoms[global_residue_index + sidechain_offset + 2] = current_atom;
                    pdb_atom_id_to_system_atom_index[atom->id] = global_residue_index + sidechain_offset + 2;
                } else if (smogTypeMap_.find(backbone_key) != smogTypeMap_.end()) {
                    int smogType = smogTypeMap_[backbone_key].first;
                    float smogRadius = smogTypeMap_[backbone_key].second;
                    mcpu_atom current_atom = mcpu_atom(chain_id, global_residue_id, n_atoms_, atom.get(), residue.get(), smogType, smogRadius);
                    if (atom->name == "N") {
                        current_atom.col_index = global_residue_index;
                        mcpu_atoms[global_residue_index] = current_atom;
                        pdb_atom_id_to_system_atom_index[atom->id] = global_residue_index;
                    } else if (atom->name == "CA") {
                        current_atom.col_index = global_residue_index + 1;
                        mcpu_atoms[global_residue_index + 1] = current_atom;
                        pdb_atom_id_to_system_atom_index[atom->id] = global_residue_index + 1;
                    } else if (atom->name == "C") {
                        current_atom.col_index = global_residue_index + sidechain_atom_count.at(residue->name) + 2;
                        mcpu_atoms[global_residue_index + sidechain_atom_count.at(residue->name) + 2] = current_atom;
                        pdb_atom_id_to_system_atom_index[atom->id] = global_residue_index + sidechain_atom_count.at(residue->name) + 2;
                    } else if (atom->name == "O") {
                        current_atom.col_index = global_residue_index + sidechain_atom_count.at(residue->name) + 3;
                        mcpu_atoms[global_residue_index + sidechain_atom_count.at(residue->name) + 3] = current_atom;
                        pdb_atom_id_to_system_atom_index[atom->id] = global_residue_index + sidechain_atom_count.at(residue->name) + 3;
                    } else if (atom->name == "OXT") {
                        // OXT is optional, but if present, it should be included                        
                        current_atom.col_index = global_residue_index + sidechain_atom_count.at(residue->name) + 4;
                        mcpu_atoms[global_residue_index + sidechain_atom_count.at(residue->name) + 4] = current_atom;
                        pdb_atom_id_to_system_atom_index[atom->id] = global_residue_index + sidechain_atom_count.at(residue->name) + 4;
                    } else if (atom->name == "OCT") {
                        // OCT is optional, but if present, it should be included
                        // OXT and OCT should not be present at the same time
                        current_atom.col_index = global_residue_index + sidechain_atom_count.at(residue->name) + 4;
                        mcpu_atoms[global_residue_index + sidechain_atom_count.at(residue->name) + 4] = current_atom;
                        pdb_atom_id_to_system_atom_index[atom->id] = global_residue_index + sidechain_atom_count.at(residue->name) + 4;
                    }
                }  else {
                    std::cerr << "Error: Smog Type not found for atom " << atom->id << std::endl;
                    continue; // Skip this atom if smog type is not found
                }
            }
        }
    }
}

void SystemMCPU::init_backbone_atoms() {
    backbone_atoms.resize(n_chains_);
    for (size_t i = 0; i < n_chains_; ++i) {
        backbone_atoms[i].resize(n_residues_per_chain[i] * 4);
    }
    for (auto& mcpu_atom: mcpu_atoms) {
        int chain_id = mcpu_atom.chain_id;
        int residue_id = mcpu_atom.residue->id; // should get local id (i.e. residue id within the chain)
        if (mcpu_atom.atom->name == "N") {
            backbone_atoms[chain_id][residue_id * 4] = mcpu_atom.col_index;
        } else if (mcpu_atom.atom->name == "CA") {
            backbone_atoms[chain_id][residue_id * 4 + 1] = mcpu_atom.col_index;
        } else if (mcpu_atom.atom->name == "C") {
            backbone_atoms[chain_id][residue_id * 4 + 2] = mcpu_atom.col_index;
        } else if (mcpu_atom.atom->name == "O") {
            backbone_atoms[chain_id][residue_id * 4 + 3] = mcpu_atom.col_index;
        }
    }
}

void SystemMCPU::init_sidechain_dihedrals() {
    // Initialize sidechain dihedrals
    // As sidechain atoms are correctly ordered in the mcpu_atoms vector
    // knowing N CA atom indices and number of dihedral angles per residue is sufficient
    number_of_sidechain_dihedrals_per_chain.resize(n_chains_);
    for (int i = 0; i < n_chains_; ++i) {
        auto& chain = protein_structure_.chains[i];
        number_of_sidechain_dihedrals_per_chain[i].resize(n_residues_per_chain[i]);
        for (auto& residue: chain->residues) {
            const std::vector<TorsionData>& torsion_data = torsion_map.at(residue->name);
            number_of_sidechain_dihedrals_per_chain[i][residue->id] = torsion_data.size();
        }
    }
}

void SystemMCPU::init_aromatic_interactions() {
    // Initialize aromatic atoms IDs
    aromatic_interactions.clear();
    for (auto& chain: protein_structure_.chains) {
        for (auto& residue: chain->residues) {
            if (aromatic_map.find(residue->name) != aromatic_map.end()) {
                auto& aromatic_data = aromatic_map.at(residue->name);
                std::array<int, 3> aromatic_atom_ids;
                for (int i = 0; i < 3; ++i) {
                    for (auto& atom: residue->atoms) {
                        if (atom->name == aromatic_data.aromatic_ring_atoms[i]) {
                            aromatic_atom_ids[i] = pdb_atom_id_to_system_atom_index[atom->id];
                            break;
                        }
                    }
                }
                aromatic_interactions.push_back({aromatic_atom_ids});
            }
        }
    }
}

void SystemMCPU::init_hydrogen_bonds() {
    hydrogen_bonds.clear();
    hydrogen_bonds.resize(n_residues_);
    int residue_id_offset = 0;
    for (auto& chain: protein_structure_.chains) {
        for (int i = 0; i < n_residues_per_chain[chain->id]; ++i) {
            auto& residue = chain->residues[i];
            auto& hbond_data = hydrogen_bonds[i + residue_id_offset];
            for (auto& atom: residue->atoms) {
                int system_atom_index = pdb_atom_id_to_system_atom_index[atom->id];
                if (atom->name == "N") {
                    hbond_data.donor_n_triangle_atoms[1] = system_atom_index;
                } else if (atom->name == "CA") {
                    hbond_data.donor_n_triangle_atoms[2] = system_atom_index;
                    hbond_data.acceptor_o_triangle_atoms[0] = system_atom_index;
                } else if (atom->name == "C") {
                    hbond_data.acceptor_o_triangle_atoms[1] = system_atom_index;
                }
            }
            if (i > 0) {
                auto& residue = chain->residues[i-1];
                for (auto& atom: residue->atoms) {
                    int system_atom_index = pdb_atom_id_to_system_atom_index[atom->id];
                    if (atom->name == "C") {
                        hbond_data.donor_n_triangle_atoms[0] = system_atom_index;
                        break;
                    }
                }
            }
            if (i < n_residues_per_chain[chain->id] - 1) {
                auto& residue = chain->residues[i+1];
                for (auto& atom: residue->atoms) {
                    int system_atom_index = pdb_atom_id_to_system_atom_index[atom->id];
                    if (atom->name == "N") {
                        hbond_data.acceptor_o_triangle_atoms[2] = system_atom_index;
                        break;
                    }
                }
            }
        }
        residue_id_offset += n_residues_per_chain[chain->id];
    }
}

void SystemMCPU::init_system_forces() {
    // Initialize forces for the system
    init_mu_potential_coefficients();
    init_torsion_parameters();
    init_aromatic_parameters();
    init_hydrogen_bond_parameters();
}

void SystemMCPU::init_mu_potential_coefficients() {
    // Initialize the coefficients for the mu potential
    // Compute mu potential parameters implementation
    Eigen::MatrixXf muPotential = forcefield_.muPotential;
    mu_potential_coefficients = Eigen::MatrixXf::Zero(n_atoms_, n_atoms_);
    for (int i = 0; i < n_atoms_; ++i) {
        int smog_type_1 = mcpu_atoms[i].smog_type;
        for (int j = i + 1; j < n_atoms_; ++j) {
            int smog_type_2 = mcpu_atoms[j].smog_type;
            mu_potential_coefficients(i, j) = muPotential(smog_type_1, smog_type_2);
        }
    }
}

void SystemMCPU::init_torsion_parameters() {
    // Compute torsional parameters implementation
    BackboneTorsionParams bbTorsions = forcefield_.bbTorsions;
    SidechainTorsionParams scTorsions = forcefield_.scTorsions;
    backbone_torsion_parameters_per_chain.resize(n_chains_);
    sidechain_torsion_parameters_per_chain.resize(n_chains_);
    residue_types_per_chain.resize(n_chains_);
    for (const auto& chain : protein_structure_.chains) {
        residue_types_per_chain[chain->id].resize(n_residues_per_chain[chain->id]);
        for (const auto& residue : chain->residues) {
            std::string resName = residue->name;
            int residue_type = residue_type_map.at(resName);
            residue_types_per_chain[chain->id][residue->id] = residue_type;
        }
    }
    for (int i = 0; i < n_chains_; ++i) {
        std::vector<int> residue_types = residue_types_per_chain[i];
        backbone_torsion_parameters_per_chain[i].resize(n_residues_per_chain[i] - 2);
        sidechain_torsion_parameters_per_chain[i].resize(n_residues_per_chain[i] - 2);
        for (int j = 1; j < n_residues_per_chain[i] - 1; ++j) {
            BackboneTorsionPerTriplet bb_torsion_params = bbTorsions[residue_types[j - 1]][residue_types[j]][residue_types[j + 1]];
            backbone_torsion_parameters_per_chain[i][j - 1] = bb_torsion_params;
            SidechainTorsionPerTriplet sc_torsion_params = scTorsions[residue_types[j - 1]][residue_types[j]][residue_types[j + 1]];
            sidechain_torsion_parameters_per_chain[i][j - 1] = sc_torsion_params;
        }
    }
}

void SystemMCPU::init_aromatic_parameters() {
    aromatic_parameters = forcefield_.aromaticPotential;
}

void SystemMCPU::init_hydrogen_bond_parameters() {
    hydrogen_bond_parameters = forcefield_.hbondPotential;
    sequence_based_hydorgen_bond_parameters = forcefield_.sequenceDependentHBondPotential;
}

int SystemMCPU::get_global_residue_index(int chain_id, int residue_id) {
    // Get the global index of a residue based on chain ID and residue ID
    if (chain_id < 0 || chain_id >= n_chains_) {
        std::cerr << "Error: Invalid chain ID " << chain_id << std::endl;
        return -1; // Invalid chain ID
    }
    int index = 0;
    for (int i = 0; i < chain_id; ++i) {
        index += n_atoms_per_chain[i]; // Sum up atoms in previous chains
    }
    for (int i = 0; i < residue_id; ++i) {
        index += protein_structure_.chains[chain_id]->residues[i]->atoms.size(); // Sum up atoms in previous residues of the same chain
    }
    return index; // Return the global index for the given chain ID
}

int SystemMCPU::get_global_residue_id(int chain_id, int residue_id) {
    // Get the global index of a residue based on chain ID and residue ID
    if (chain_id < 0 || chain_id >= n_chains_) {
        std::cerr << "Error: Invalid chain ID " << chain_id << std::endl;
        return -1; // Invalid chain ID
    }
    int id = 0;
    for (int i = 0; i < chain_id; ++i) {
        id += n_residues_per_chain[i]; // Sum up atoms in previous chains
    }
    id += residue_id; // Add the residue ID in the current chain
    return id; // Return the global residue id
}