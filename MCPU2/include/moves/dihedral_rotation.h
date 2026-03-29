#ifndef GLOBAL_ROTATION_H
#define GLOBAL_ROTATION_H

#include <Eigen/Dense>

class SystemMCPU;

class DihedralRotation {
   public:
    DihedralRotation() = default;
    // Rotate dihedral angle
    static Eigen::Matrix3Xf rotate_dihedral(Eigen::Matrix3Xf       &positions,
                                            const std::vector<int> &atom_ids, float angle);
    // Rotate phi angle
    static Eigen::Matrix3Xf rotate_phi(Eigen::Matrix3Xf &positions, int &chain_id, int &residue_id);
    // Rotate psi angle
    static Eigen::Matrix3Xf rotate_psi(Eigen::Matrix3Xf &positions, int &chain_id, int &residue_id);
    // Rotate chi angle
    static Eigen::Matrix3Xf rotate_chi(Eigen::Matrix3Xf &positions, int &chain_id, int &residue_id);
    // Generate residue_id to atom_ids map
    void generate_residue_id_to_rotating_atom_ids_map(const std::vector<int> &backbone_atom_ids);

    // Rotation Atom IDs
    static std::vector<std::unordered_map<int, std::vector<int>>> phi_rotation_atom_ids;
    static std::vector<std::unordered_map<int, std::vector<int>>> psi_rotation_atom_ids;
    static std::vector<std::unordered_map<int, std::vector<int>>> chi_rotation_atom_ids;

   private:
    SystemMCPU &system_;
};

#endif  // GLOBAL_ROTATION_H