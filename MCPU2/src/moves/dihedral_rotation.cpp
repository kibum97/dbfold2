#include "moves/dihedral_rotation.h"
#include "utils/rotation_utils.h"
#include "utils/random_utils.h"
#include "system_mcpu.h"
#include <numeric>

Eigen::Matrix3Xf DihedralRotation::rotate_dihedral(Eigen::Matrix3Xf& positions, const std::vector<int>& atom_ids, float angle) {
    // Set origin
    Eigen::Vector3f origin = positions.col(atom_ids[0]);
    // Set axis
    Eigen::Vector3f axis = positions.col(atom_ids[1]) - positions.col(atom_ids[0]);
    axis.normalize();
    // Set rotation
    Eigen::AngleAxisf rotation = Eigen::AngleAxisf(angle, axis);
    for (auto atomid = atom_ids.begin() + 2; atomid != atom_ids.end(); ++atomid) {
        size_t atom_id = *atomid;
        positions.col(atom_id) = rotate_vector(positions.col(atom_id), origin, rotation);
    }
    return positions;
}

Eigen::Matrix3Xf DihedralRotation::rotate_phi(Eigen::Matrix3Xf& positions, int& chain_id, int& residue_id) {
    float angle = generate_gaussian(0, 1);
    return rotate_dihedral(positions, phi_rotation_atom_ids[chain_id][residue_id], angle);
}

Eigen::Matrix3Xf DihedralRotation::rotate_psi(Eigen::Matrix3Xf& positions, int& chain_id, int& residue_id) {
    float angle = generate_gaussian(0, 1);
    return rotate_dihedral(positions, psi_rotation_atom_ids[chain_id][residue_id], angle);
}

Eigen::Matrix3Xf DihedralRotation::rotate_chi(Eigen::Matrix3Xf& positions, int& chain_id, int& residue_id) {
    float angle = generate_gaussian(0, 1);
    return rotate_dihedral(positions, chi_rotation_atom_ids[chain_id][residue_id], angle);
}

/*
Eigen::Matrix3Xf DihedralRotation::sidechain_rotation(size_t resID, Eigen::Matrix3Xd& positions, std::unordered_map<size_t, std::vector<RotamerData>> rotamerIDMap, std::unordered_map<size_t, std::vector<TorsionIDData>> rotatingSCAtomsMap, std::vector<std::array<double, 4>> sidechain_torsions) {
    // Get rotation angles from rotamer data
    std::vector<size_t> rotamerIDs = {};
    std::vector<double> weights = {};
    size_t index = 0;
    for (const auto& rotamerData : rotamerIDMap[resID]) {
        rotamerIDs.push_back(index);
        weights.push_back(rotamerData.probability);
        ++index;
    }
    if (rotamerIDs.size() > 0) {
        RotamerData selected_rotamer = rotamerIDMap[resID][selectElementWithWeights(rotamerIDs, weights)];
        // Perform sidechain rotation
        for (const auto& torsionIDData : rotatingSCAtomsMap[resID]) {
            int torsion_id = torsionIDData.torsion_id;
            std::vector<size_t> atomIDs;
            atomIDs.push_back(torsionIDData.torsion_atomIDs[1]);
            atomIDs.push_back(torsionIDData.torsion_atomIDs[2]);
            for (size_t atomID : torsionIDData.rotating_atomIDs) {
                atomIDs.push_back(atomID);
            }
            Eigen::Vector3d axis = positions.col(atomIDs[1]) - positions.col(atomIDs[0]);
            axis.normalize();
            double angle = generateGaussian(selected_rotamer.chiMeans[torsion_id],selected_rotamer.chiStdDevs[torsion_id]) - sidechain_torsions[resID][torsion_id];
            positions = partialRotation(positions, atomIDs, angle, axis);
        }
    }
    return positions;
}
*/

void DihedralRotation::generate_residue_id_to_rotating_atom_ids_map(const std::vector<int>& backbone_atom_ids) {
    // Backbone atom indices of given chain
    for (int chain_id = 0; chain_id < system_.n_residues_per_chain.size(); ++chain_id) {
        std::vector<int> backbone_atom_ids = system_.backbone_atoms[chain_id];
        int n_residues = system_.n_residues_per_chain[chain_id];
        int n_atoms = system_.n_atoms_per_chain[chain_id];
        for (int i = 0; i < n_residues; ++i) {
            int n_index = backbone_atom_ids[i * 4];
            int ca_index = backbone_atom_ids[i * 4 + 1];
            int c_index = backbone_atom_ids[i * 4 + 2];
            for (int phi_index = n_index; phi_index < n_atoms; ++phi_index) {
                phi_rotation_atom_ids[chain_id][i].push_back(phi_index);
            }
            for (int psi_index = ca_index; psi_index < n_atoms; ++psi_index) {
                psi_rotation_atom_ids[chain_id][i].push_back(psi_index);
            }
            for (int chi_index = n_index; chi_index < c_index; ++chi_index) {
                // TODO: chi angle of later chis
                chi_rotation_atom_ids[chain_id][i].push_back(chi_index);
            }
        }
    }
}