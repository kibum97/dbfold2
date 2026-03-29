#ifndef HBOND_H
#define HBOND_H

#include <Eigen/Dense>

/*

bool check_hydrogen_bond(size_t donorResID, std::array<size_t, 7> donorAtomIDs,
                         size_t acceptorResID, std::array<size_t, 8> acceptorAtomIDs,
                         Eigen::Matrix3Xf positions);
std::array<float, 3> compute_distance_features(size_t donorResID, std::array<size_t, 7>
donorAtomIDs, size_t acceptorResID, std::array<size_t, 8> acceptorAtomIDs, Eigen::Matrix3Xf
positions); std::array<float, 4> compute_dihedral_features(size_t donorResID, std::array<size_t, 7>
donorAtomIDs, size_t acceptorResID, std::array<size_t, 8> acceptorAtomIDs, Eigen::Matrix3Xf
positions); Eigen::Vector3f infer_hydrogen_position(const Eigen::Vector3f& c, const Eigen::Vector3f&
n, const Eigen::Vector3f& ca); int compute_hbond_index_from_angles(int ss, std::array<float, 3>
bisect_angles, std::array<float, 3> plane_angles);
*/
std::array<float, 3> compute_distance_features(int donor_residue_index, int acceptor_residue_index,
                                               Eigen::MatrixXf ca_ca_distances);
std::array<int, 3>   compute_angular_features(
      int donor_chain_id, int donor_residue_id, int acceptor_chain_id, int acceptor_residue_id,
      std::vector<std::vector<Eigen::Vector3f>> &backbone_vectors);
std::array<float, 4> compute_dihedral_features(
    int donor_chain_id, int donor_residue_id, int acceptor_chain_id, int acceptor_residue_id,
    std::vector<std::vector<float>> &backbone_phi_angles,
    std::vector<std::vector<float>> &backbone_psi_angles);
int get_residue_difference(int donor_chain_id, int donor_residue_id, int acceptor_chain_id,
                           int acceptor_residue_id);

bool check_hydrogen_bond(int donor_chain_id, int donor_residue_id, int donor_global_residue_index,
                         int acceptor_chain_id, int acceptor_residue_id,
                         int                              acceptor_global_residue_index,
                         std::vector<std::vector<float>> &backbone_phi_angles,
                         std::vector<std::vector<float>> &backbone_psi_angles,
                         Eigen::MatrixXf                  ca_ca_distances);
int  compute_secondary_structure(Eigen::Vector3f donor_ca_vector,
                                 Eigen::Vector3f acceptor_ca_vector);
int  compute_hbond_index(int ss, std::array<int, 3> bisect_indices,
                         std::array<int, 3> plane_indices);

#endif  // HBOND_H