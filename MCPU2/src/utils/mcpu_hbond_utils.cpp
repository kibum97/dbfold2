#include "utils/mcpu_hbond_utils.h"
#include "utils/geometry_utils.h"
#include <iostream>
#include <algorithm>
#include <cstdlib>

/*
Eigen::Vector3f infer_hydrogen_position(const Eigen::Vector3f& c, const Eigen::Vector3f& n, const Eigen::Vector3f& ca) {
    // Infer the hydrogen position
    Eigen::Vector3f b1 = ca - n;
    Eigen::Vector3f b2 = c - n;
    Eigen::Vector3f h_norm = (b1 + b2).normalized();
    h_norm *= -1.0;
    Eigen::Vector3f hydrogen = n + h_norm;
    return hydrogen;
}

int compute_hbond_index_from_angles(int ss, std::array<float, 3> bisect_angles, std::array<float, 3> plane_angles) {
    // Compute hydrogen bond index from hbond state
    std::array<int, 3> bisect_indices;
    std::array<int, 3> plane_indices;
    for (size_t i = 0; i < 3; ++i) {
        bisect_indices[i] = static_cast<int>(bisect_angles[i]  * (180.0 / M_PI) / 20.);
        plane_indices[i] = static_cast<int>(plane_angles[i] * (180.0 / M_PI) / 20.);
    }
    // Compute hydrogen bond index
    return bisect_indices[2] 
       + plane_indices[2] * 9 
       + bisect_indices[1]  * 9 * 9 
       + plane_indices[1] * 9 * 9 * 9 
       + bisect_indices[0]  * 9 * 9 * 9 * 9 
       + plane_indices[0] * 9 * 9 * 9 * 9 * 9 
       + ss * 9 * 9 * 9 * 9 * 9 * 9;
}
*/



std::array<float, 3> compute_distance_features(int donor_residue_index, int acceptor_residue_index, Eigen::MatrixXf ca_ca_distances) {
    // Secondary structure specific - CA distances
    float ca_dist_21 = ca_ca_distances(donor_residue_index - 1, acceptor_residue_index);
    float ca_dist_10 = ca_ca_distances(donor_residue_index, acceptor_residue_index + 1);
    float ca_dist_20 = ca_ca_distances(donor_residue_index -1, acceptor_residue_index + 1);
    float ca_dist_11 = ca_ca_distances(donor_residue_index, acceptor_residue_index);
    // Compute minimums
    float min1 = std::min(ca_dist_21, ca_dist_20);
    float min2 = std::min(ca_dist_10, ca_dist_11);
    float min3 = std::min(min1, min2);
    return {min1, min2, min3};
}

std::array<int, 3> compute_angular_features(int donor_chain_id, int donor_residue_id,
                                            int acceptor_chain_id, int acceptor_residue_id,
                                            std::vector<std::vector<Eigen::Vector3f>>& backbone_vectors) {
    float angle1 = compute_angle(backbone_vectors[donor_chain_id][donor_residue_id * 3],
                                        backbone_vectors[acceptor_chain_id][acceptor_residue_id * 3]);
    float angle2 = compute_angle(backbone_vectors[donor_chain_id][(donor_residue_id - 1) * 3],
                                        backbone_vectors[acceptor_chain_id][(acceptor_residue_id + 1) * 3]);
    float angle3 = compute_angle(backbone_vectors[donor_chain_id][(donor_residue_id - 1) * 3 + 2],
                                        backbone_vectors[acceptor_chain_id][acceptor_residue_id * 3 + 1]);
    int angle1_index = static_cast<int>(angle1 * (180.0 / M_PI) / 20.);
    int angle2_index = static_cast<int>(angle2 * (180.0 / M_PI) / 20.);
    int angle3_index = static_cast<int>(angle3 * (180.0 / M_PI) / 20.);
    return {angle1_index, angle2_index, angle3_index};
}

std::array<float, 4> compute_dihedral_features(int donor_chain_id, int donor_residue_id,
                                               int acceptor_chain_id, int acceptor_residue_id,
                                               std::vector<std::vector<float>>& backbone_phi_angles,
                                               std::vector<std::vector<float>>& backbone_psi_angles) {
    // Secondary structure specific - dihedrals
    float donor_phi = backbone_phi_angles[donor_chain_id][donor_residue_id - 1];
    float donor_psi = backbone_psi_angles[donor_chain_id][donor_residue_id - 2];   
    float acceptor_phi = backbone_phi_angles[acceptor_chain_id][acceptor_residue_id];
    float acceptor_psi = backbone_psi_angles[acceptor_chain_id][acceptor_residue_id - 1];
    // Adjust range
    donor_phi += 180.0;
    donor_psi += 180.0;
    acceptor_phi += 180.0;
    acceptor_psi += 180.0;
    return {donor_phi, donor_psi, acceptor_phi, acceptor_psi};
}

bool check_hydrogen_bond(int donor_chain_id, int donor_residue_id, int donor_global_residue_index,
                         int acceptor_chain_id, int acceptor_residue_id, int acceptor_global_residue_index,
                         std::vector<std::vector<float>>& backbone_phi_angles,
                         std::vector<std::vector<float>>& backbone_psi_angles,
                         Eigen::MatrixXf ca_ca_distances) {
    int residue_diff = get_residue_difference(donor_chain_id, donor_residue_id, acceptor_chain_id, acceptor_residue_id);
    // Check if a hydrogen bond is formed
    if (residue_diff < 4) {
        return false;
    } else if (residue_diff == 4) {
        std::array<float, 3> distance_features = compute_distance_features(donor_global_residue_index, acceptor_global_residue_index, ca_ca_distances);
        if ((distance_features[0] > 5.8 * 5.8) || (distance_features[1] > 5.8 * 5.8) || (distance_features[2] > 5.5 * 5.5)) {
            return false;
        }
    } else {
        std::array<float, 3> distance_features = compute_distance_features(donor_global_residue_index, acceptor_global_residue_index, ca_ca_distances);
        if ((distance_features[0] > 6.0 * 6.0) || (distance_features[1] > 6.0 * 6.0) || (distance_features[2] > 5.4 * 5.4)) {
            return false;
        }
        std::array<float, 4> dihedral_features = compute_dihedral_features(donor_chain_id, donor_residue_id, acceptor_chain_id, acceptor_residue_id,
                                                                           backbone_phi_angles, backbone_psi_angles);
        if ((dihedral_features[0] > 150.) || (dihedral_features[2] > 150.)) {
            return false;
        }
        if ((dihedral_features[1] > 30.) && (dihedral_features[1] < 210.) || (dihedral_features[3] > 30.) && (dihedral_features[3] < 210.)) {
            return false;
        }
    }
    return true;
}

int compute_secondary_structure(Eigen::Vector3f donor_ca_vector, Eigen::Vector3f acceptor_ca_vector) {
    // If residue number difference between donor and acceptor is larger than 4
    // call this function to compute secondary structure label
    int secondary_structure_label = 0;
    float angle = compute_angle(donor_ca_vector, acceptor_ca_vector);
    angle = angle * 180.0 / M_PI;
    if (angle < 90.0) {
        secondary_structure_label = 1;
    } else {
        secondary_structure_label = 2;
    }
    return secondary_structure_label;
}

int get_residue_difference(int donor_chain_id, int donor_residue_id,
                           int acceptor_chain_id, int acceptor_residue_id) {
    // Get the residue difference
    if (donor_chain_id == acceptor_chain_id) {
        return std::abs(donor_residue_id - acceptor_residue_id);
    } else {
        return 100; // Different chains
    }
}

int compute_hbond_index(int ss, std::array<int, 3> bisect_indices, std::array<int, 3> plane_indices) {  
    // Compute hydrogen bond index
    return bisect_indices[2] 
       + plane_indices[2] * 9 
       + bisect_indices[1]  * 9 * 9 
       + plane_indices[1] * 9 * 9 * 9 
       + bisect_indices[0]  * 9 * 9 * 9 * 9 
       + plane_indices[0] * 9 * 9 * 9 * 9 * 9 
       + ss * 9 * 9 * 9 * 9 * 9 * 9;
}