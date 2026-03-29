#include "system_mcpu.h"
#include "moves/concerted_rotation.h"
#include "tripeptide_closure/tripep_closure.h"
#include "utils/geometry_utils.h"
#include "utils/random_utils.h"

#define STATUS stdout

ConcertedRotation::ConcertedRotation(SystemMCPU& system) :
    system_(system) {
};

Eigen::Matrix3Xf ConcertedRotation::rotate(Eigen::Matrix3Xf& positions, int residue_id, float angle) {
    // Perform concerted rotation
    int n_solutions = 0;
    std::vector<size_t> involved_atom_ids;
    std::array<Eigen::Vector3d, 15> old_backbone_positions;
    std::array<Eigen::Vector3d, 15> new_backbone_positions;
    std::vector<double> delta_dihedrals;
    Eigen::Matrix3Xd new_positions = positions;
    double r_n[3][3], r_a[3][3], r_c[3][3];
    double r_soln_n[max_soln][3][3], r_soln_a[max_soln][3][3], r_soln_c[max_soln][3][3];
    double b_len[6], b_ang[7], t_ang[2];
    // Get the list of involved atoms
    involved_atom_ids = rotation_atom_ids[chain_id][residue_id];
    for (size_t i = 0; i < involved_atom_ids.size(); ++i) {
        old_backbone_positions[i] = positions.col(involved_atom_ids[i]);
        new_backbone_positions[i] = positions.col(involved_atom_ids[i]);
    }
    // Store the positions of the atoms into correct format
    for (size_t i = 0; i < resIDs.size(); ++i) {
        for (int k = 0; k < 3; ++k) {
                r_n[i][k] = positions(k, backbo);
                r_a[i][k] = positions(k, atomid);
                r_c[i][k] = positions(k, atomid);
        }
    }
    // Initialize loop closure
    for (int bi = 0; bi < 6; ++bi) {
        b_len[bi] = compute_length(positions.col(involved_atom_ids[bi+4]), positions.col(involved_atom_ids[bi+5]));
    }
    for (int ai = 0; ai < 7; ++ai) {
        b_ang[ai] = compute_bond_angle(
            positions.col(involved_atom_ids[ai+3]),
            positions.col(involved_atom_ids[ai+4]),
            positions.col(involved_atom_ids[ai+5])
        );
    }
    t_ang[0] = pi;
    t_ang[1] = pi;    
    // Solve the tripeptide closure problem
    initialize_loop_closure(b_len, b_ang, t_ang);
    solve_3pep_poly(r_n[0], r_a[0], r_a[2], r_c[2], r_soln_n, r_soln_a, r_soln_c, &n_solutions);

    // Select random index for the solution
    if (n_solutions <= 0) {
        return std::make_tuple(positions, n_solutions);
    } else {
        int random_number = getRandomNumber(n_solutions);
        random_number = 0;
        
        // Store the positions of the atoms back into the matrix
        for (size_t i = 0; i < resIDs.size(); ++i) {
            for (int l = 0; l < 3; ++l) {
                size_t atomid = backboneAtomIDs[resIDs[i]][l];
                Eigen::Vector3d solution;
                if (l == 0) {
                    solution = Eigen::Vector3d(r_soln_n[random_number][i][0], r_soln_n[random_number][i][1], r_soln_n[random_number][i][2]);
                } else if (l == 1) {
                    solution = Eigen::Vector3d(r_soln_a[random_number][i][0], r_soln_a[random_number][i][1], r_soln_a[random_number][i][2]);
                } else if (l == 2) {
                    solution = Eigen::Vector3d(r_soln_c[random_number][i][0], r_soln_c[random_number][i][1], r_soln_c[random_number][i][2]);
                }
                new_backbone_positions[i*3 + l + 3] = solution;
            }
        }

        for (size_t i = 3; i < involved_atom_ids.size() - 4; ++i) {
            if (i % 3 == 0) {
                // Rotate sidechain R
                size_t resID = involved_res_ids[i/3];
                double delta_dihedral = compute_dihedral_angle(
                    new_positions.col(involved_atom_ids[i-1]),
                    new_positions.col(involved_atom_ids[i]),
                    new_positions.col(involved_atom_ids[i+1]),
                    new_backbone_positions[i+2]
                ) - compute_dihedral_angle(
                    new_positions.col(involved_atom_ids[i-1]),
                    new_positions.col(involved_atom_ids[i]),
                    new_positions.col(involved_atom_ids[i+1]),
                    new_positions.col(involved_atom_ids[i+2])
                );
                std::vector<size_t> rotating_IDs = {involved_atom_ids[i], involved_atom_ids[i+1]};
                rotating_IDs.insert(rotating_IDs.end(), involved_atom_ids.begin() + i + 2, involved_atom_ids.end());
                for (size_t rid = resID; rid < resIDs[2] + 1; ++rid) {
                    rotating_IDs.insert(rotating_IDs.end(), sidechainAtomIDs[rid].begin(), sidechainAtomIDs[rid].end());
                    rotating_IDs.insert(rotating_IDs.end(), backboneAtomIDs[rid][3]);
                }
                Eigen::Vector3d axis = new_positions.col(rotating_IDs[1]) - new_positions.col(rotating_IDs[0]);
                axis.normalize();
                new_positions = partialRotation(new_positions, rotating_IDs, delta_dihedral, axis);

            } else if (i % 3 == 1) {
                // Rotate backbone O
                size_t resID = involved_res_ids[i/3];
                double delta_dihedral = compute_dihedral_angle(
                    new_positions.col(involved_atom_ids[i-1]),
                    new_positions.col(involved_atom_ids[i]),
                    new_positions.col(involved_atom_ids[i+1]),
                    new_backbone_positions[i+2]
                ) - compute_dihedral_angle(
                    new_positions.col(involved_atom_ids[i-1]),
                    new_positions.col(involved_atom_ids[i]),
                    new_positions.col(involved_atom_ids[i+1]),
                    new_positions.col(involved_atom_ids[i+2])
                );
                std::vector<size_t> rotating_IDs = {involved_atom_ids[i] , involved_atom_ids[i+1], backboneAtomIDs[resID][3]};
                rotating_IDs.insert(rotating_IDs.end(), involved_atom_ids.begin() + i + 2, involved_atom_ids.end());
                for (size_t rid = resID + 1; rid < resIDs[2] + 1; ++rid) {
                    rotating_IDs.insert(rotating_IDs.end(), sidechainAtomIDs[rid].begin(), sidechainAtomIDs[rid].end());
                    rotating_IDs.insert(rotating_IDs.end(), backboneAtomIDs[rid][3]);
                }
                Eigen::Vector3d axis = new_positions.col(rotating_IDs[1]) - new_positions.col(rotating_IDs[0]);
                axis.normalize();
                new_positions = partialRotation(new_positions, rotating_IDs, delta_dihedral, axis);
            } else {
                size_t resID = involved_res_ids[i/3];
                double delta_dihedral = compute_dihedral_angle(
                    new_positions.col(involved_atom_ids[i-1]),
                    new_positions.col(involved_atom_ids[i]),
                    new_positions.col(involved_atom_ids[i+1]),
                    new_backbone_positions[i+2]
                ) - compute_dihedral_angle(
                    new_positions.col(involved_atom_ids[i-1]),
                    new_positions.col(involved_atom_ids[i]),
                    new_positions.col(involved_atom_ids[i+1]),
                    new_positions.col(involved_atom_ids[i+2])
                );
                std::vector<size_t> rotating_IDs = {involved_atom_ids[i] , involved_atom_ids[i+1]};
                rotating_IDs.insert(rotating_IDs.end(), involved_atom_ids.begin() + i + 2, involved_atom_ids.end());
                for (size_t rid = resID + 1; rid < resIDs[2] + 1; ++rid) {
                    rotating_IDs.insert(rotating_IDs.end(), sidechainAtomIDs[rid].begin(), sidechainAtomIDs[rid].end());
                    rotating_IDs.insert(rotating_IDs.end(), backboneAtomIDs[rid][3]);
                }
                Eigen::Vector3d axis = new_positions.col(rotating_IDs[1]) - new_positions.col(rotating_IDs[0]);
                axis.normalize();
                new_positions = partialRotation(new_positions, rotating_IDs, delta_dihedral, axis);
            }
        }
    }
    return new_positions;
}

void ConcertedRotation::generate_pivot_atom_id_to_rotating_atom_ids_map(const std::vector<int>& backbone_atom_ids) {
    // Generate the rotating atom map based on the system
    for (int chain_id = 0; chain_id < system_.n_residues_per_chain.size(); ++chain_id) {
        for (int res_id = 0; res_id < system_.n_residues_per_chain[chain_id] - 2; ++res_id) {
            int n_index = backbone_atom_ids[res_id * 4];
            for (int atom_index = n_index; atom_index < backbone_atom_ids[(res_id + 3) * 4]; ++atom_index) {
                rotation_atom_ids[chain_id][res_id].push_back(atom_index);
            }
        }
    }
}