#include "utils/hbond.h"

int computeSecondaryStructure(size_t donorResID, std::array<size_t, 7> donorAtomIDs,
                              size_t acceptorResID, std::array<size_t, 8> acceptorAtomIDs,
                              Eigen::Matrix3Xd positions) {
    size_t residue_diff = std::abs(static_cast<int>(acceptorResID - donorResID));
    int secondary_structure_label = 0;
    if (residue_diff == 4) {
        secondary_structure_label = 0;
    } else if (residue_diff > 4) {
        double angle = computeAngle(positions.col(donorAtomIDs[6]) - positions.col(donorAtomIDs[4]),
                                    positions.col(acceptorAtomIDs[5]) - positions.col(acceptorAtomIDs[7]));
        angle = angle * 180.0 / M_PI;
        if (angle < 90.0) {
            secondary_structure_label = 1;
        } else {
            secondary_structure_label = 2;
        }
    }
    return secondary_structure_label;
}

bool checkHydrogenBond(size_t donorResID, std::array<size_t, 7> donorAtomIDs,
                       size_t acceptorResID, std::array<size_t, 8> acceptorAtomIDs,
                       Eigen::Matrix3Xd positions) {
    size_t residue_diff = std::abs(static_cast<int>(acceptorResID - donorResID));
    // Check if a hydrogen bond is formed
    Eigen::Vector3d donor_hydrogen = inferHydrogenPosition(positions.col(donorAtomIDs[1]), positions.col(donorAtomIDs[0]), positions.col(donorAtomIDs[2]));
    double donor_acceptor_dist = (donor_hydrogen - positions.col(acceptorAtomIDs[0])).norm();
    if (donor_acceptor_dist > 2.5) {
        return false;
    } else {
        std::array<double, 3> distance_features = computeDistanceFeatures(donorResID, donorAtomIDs, acceptorResID, acceptorAtomIDs, positions);
        if (residue_diff == 4) {
            if ((distance_features[0] > 5.8) || (distance_features[1] > 5.8) || (distance_features[2] > 5.5)) {
                return false;
            }
        } else if (residue_diff > 4) {
            if ((distance_features[0] > 6.0) || (distance_features[1] > 6.0) || (distance_features[2] > 5.4)) {
                return false;
            }
        }
        std::array<double, 4> dihedral_features = computeDihedralFeatures(donorResID, donorAtomIDs, acceptorResID, acceptorAtomIDs, positions);
        if (residue_diff > 4) {
            // don't form long-range h-bonds if we occupy quadrants on the right of ramachandran plot
            if ((dihedral_features[0] > 150.) || (dihedral_features[2] > 150.)) {
                return false;
            }
            if ((dihedral_features[1] > 30.) && (dihedral_features[1] < 210.) || (dihedral_features[3] > 30.) && (dihedral_features[3] < 210.)) {
                return false;
            }
        }
        return true;
    }
}

std::array<double, 3> computeDistanceFeatures(size_t donorResID, std::array<size_t, 7> donorAtomIDs,
                                              size_t acceptorResID, std::array<size_t, 8> acceptorAtomIDs,
                                              Eigen::Matrix3Xd positions) {
    // Secondary structure specific - CA distances
    double ca_dist_21 = (positions.col(donorAtomIDs[4]) - positions.col(acceptorAtomIDs[2])).norm();
    double ca_dist_10 = (positions.col(donorAtomIDs[1]) - positions.col(acceptorAtomIDs[5])).norm();
    double ca_dist_20 = (positions.col(donorAtomIDs[4]) - positions.col(acceptorAtomIDs[5])).norm();
    double ca_dist_11 = (positions.col(donorAtomIDs[1]) - positions.col(acceptorAtomIDs[2])).norm();
    // Compute minimums
    double min1 = std::min(ca_dist_21, ca_dist_20);
    double min2 = std::min(ca_dist_10, ca_dist_11);
    double min3 = std::min(min1, min2);
    return {min1, min2, min3};
}

std::array<double, 4> computeDihedralFeatures(size_t donorResID, std::array<size_t, 7> donorAtomIDs,
                                              size_t acceptorResID, std::array<size_t, 8> acceptorAtomIDs,
                                              Eigen::Matrix3Xd positions) {
    // Secondary structure specific - dihedrals
    double donor_phi = computeDihedralAngle(positions.col(donorAtomIDs[2]), positions.col(donorAtomIDs[0]), positions.col(donorAtomIDs[1]), positions.col(donorAtomIDs[3]));
    double donor_psi = computeDihedralAngle(positions.col(donorAtomIDs[5]), positions.col(donorAtomIDs[4]), positions.col(donorAtomIDs[2]), positions.col(donorAtomIDs[0]));
    double acceptor_phi = computeDihedralAngle(positions.col(acceptorAtomIDs[1]), positions.col(acceptorAtomIDs[3]), positions.col(acceptorAtomIDs[5]), positions.col(acceptorAtomIDs[6]));
    double acceptor_psi = computeDihedralAngle(positions.col(acceptorAtomIDs[4]), positions.col(acceptorAtomIDs[2]), positions.col(acceptorAtomIDs[1]), positions.col(acceptorAtomIDs[3]));
    // Make angles to be in degree
    donor_phi = donor_phi * 180.0 / M_PI;
    donor_psi = donor_psi * 180.0 / M_PI;
    acceptor_phi = acceptor_phi * 180.0 / M_PI;
    acceptor_psi = acceptor_psi * 180.0 / M_PI;
    // Adjust range
    donor_phi += 180.0;
    donor_psi += 180.0;
    acceptor_phi += 180.0;
    acceptor_psi += 180.0;
    return {donor_phi, donor_psi, acceptor_phi, acceptor_psi};
}

Eigen::Vector3d inferHydrogenPosition(const Eigen::Vector3d& c, const Eigen::Vector3d& n, const Eigen::Vector3d& ca) {
    // Infer the hydrogen position
    Eigen::Vector3d b1 = ca - n;
    Eigen::Vector3d b2 = c - n;
    Eigen::Vector3d h_norm = (b1 + b2).normalized();
    h_norm *= -1.0;
    Eigen::Vector3d hydrogen = n + h_norm;
    return hydrogen;
}

int computeHBondindexfromAngles(int ss, std::array<double, 3> bisect_angles, std::array<double, 3> plane_angles) {
    // Compute hydrogen bond index from hbond state
    std::array<int, 3> bisect_indices;
    std::array<int, 3> plane_indices;
    for (size_t i = 0; i < 3; ++i) {
        bisect_indices[i] = static_cast<int>(bisect_angles[i] / 20.);
        plane_indices[i] = static_cast<int>(plane_angles[i] / 20.);
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

int computeHBondindex(int ss, std::array<int, 3> bisect_indices, std::array<int, 3> plane_indices) {  
    // Compute hydrogen bond index
    return bisect_indices[2] 
       + plane_indices[2] * 9 
       + bisect_indices[1]  * 9 * 9 
       + plane_indices[1] * 9 * 9 * 9 
       + bisect_indices[0]  * 9 * 9 * 9 * 9 
       + plane_indices[0] * 9 * 9 * 9 * 9 * 9 
       + ss * 9 * 9 * 9 * 9 * 9 * 9;
}