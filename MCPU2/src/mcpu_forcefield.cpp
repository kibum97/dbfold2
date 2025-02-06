#include "mcpu_forcefield.h"
#include "read_config.h"
#include "utils/hbond.h"
#include <fstream>
#include <iostream>

//TODO: Need to read hbond_file from config file - Need to add this var

int USE_SMOG_TYPE = 1;
extern std::string potential_file;
extern std::string atom_type_file;
extern std::string hbond_file;
extern std::string hydrogen_bonding_data;
extern std::string all_triplet_file;
extern std::string all_sctorsion_file;

MCPUForceField::MCPUForceField() {
    // Initialize force field implementation
    loadMuPotential(potential_file);
    loadSmogTypeMap(atom_type_file);
    loadHbondPotential(hbond_file);
    loadSequenceDependentHBondPotential(seq_hbond_file);
    for (int i; i < 20; i++) {
        for (int j; j < 20; j++) {
            for (int k; k < 20; k++) {
                for (int l; l < 6; l++) {
                    for (int m; m < 6; m++) {
                        for (int n; n < 6; n++) {
                            for (int o; o < 6; o++) {
                                bbTorsions[i][j][k][l][m][n][o] = 1000;
                                scTorsions[i][j][k][l][m][n][o] = 1000;
                            }
                        }
                    }
                }
            }
        }
    }
    loadTorsionalPotential(all_triplet_file, all_sctorsion_file);
}

MCPUForceField::~MCPUForceField() {
    // Destructor implementation
}

void MCPUForceField::loadMuPotential(const std::string& filename) {
    // Load mu potential implementation
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open mu potential file.");
    }

    // Find the matrix dimensions dynamically
    std::vector<std::tuple<int, int, double>> data;
    int maxRow = 0, maxCol = 0;
    int row, col;
    double value;

    while (file >> row >> col >> value) {
        data.emplace_back(row, col, value);
        if (row > maxRow) maxRow = row;
        if (col > maxCol) maxCol = col;
    }
    file.close();

    // Create Eigen matrix with the appropriate size
    Eigen::MatrixXd muParamMatrix = Eigen::MatrixXd::Zero(maxRow + 1, maxCol + 1);

    // Populate the matrix
    for (const auto& [rowIdx, colIdx, val] : data) {
        muParamMatrix(rowIdx, colIdx) = val;
        muParamMatrix(colIdx, rowIdx) = val;
    }
    muPotential = muParamMatrix;
}

void MCPUForceField::loadSmogTypeMap(const std::string& filename) {
    // Load smog type map implementation
    std::unordered_map<std::string, int> smogtypeMap;

    // Open the file
    std::ifstream inputFile(filename);
    if (!inputFile.is_open()) {
        throw std::runtime_error("Could not open SMOG type file.");
    }

    // Read the file line by line
    std::string line;
    while (std::getline(inputFile, line)) {
        std::istringstream ss(line);
        std::string atomName, residueName;
        int smogType;

        // Parse the three columns
        if (ss >> atomName >> residueName >> smogType) {
            // Concatenate the first wo columns to fdonor_orm the key
            std::string key = residueName + "_" + atomName;
            // Insert into the map
            smogtypeMap[key] = smogType;
        } else {
            std::cerr << "Error: Malformed line: " << line << std::endl;
        }
    }

    // Close the file
    inputFile.close();

    smogTypeMap = smogtypeMap;
}

void MCPUForceField::loadHbondPotential(const std::string& filename) {
    // Load hydrogen bond potential implementation
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open Hbond file.");
    }

    // Find the matrix dimensions dynamically
    int ss;
    int plane_angle_ca, bisect_angle_ca, plane_angle_ca_off, bisect_angle_ca_off;
    int bisect_angle_strand, plane_angle_strand;
    double value;
    int index;

    // Store the value in the map
    while (file >> ss >> plane_angle_ca >> bisect_angle_ca >> plane_angle_ca_off >> bisect_angle_ca_off >> plane_angle_strand >> bisect_angle_strand >> value) {    
        std::array <int, 3> bisect_indices = {bisect_angle_ca, bisect_angle_ca_off, bisect_angle_strand};
        std::array <int, 3> plane_indices = {plane_angle_ca, plane_angle_ca_off, plane_angle_strand};
        index = computeHBondindex(ss, bisect_indices, plane_indices);
        hbondPotential[index] = value;
    }
    file.close();
}

void MCPUForceField::loadSequenceDependentHBondPotential(const std::string& filename) {
    // Load sequence dependent hydrogen bond potential implementation
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open sequence dependent Hbond file.");
    }
    // Read and store values
    int secondary_structure_label, donor_resType, acceptor_resType;
    double value;
    int obs, tot_obs;
    while (file >> secondary_structure_label >> donor_resType >> acceptor_resType >> value >> obs >> tot_obs) {
        std::cout << "Value from file - HBOND SEQ DEP : " << value << std::endl;
        sequenceDependentHBondPotential[secondary_structure_label][donor_resType][acceptor_resType] = (-1.) * value;
    }
    file.close();
}

void MCPUForceField::loadTorsionalPotential(const std::string& triple_filename, const std::string& sctorsion_filename) {
    // Load backbone torsion parameters
    std::ifstream triple_file(triple_filename);
    if (!triple_file.is_open()) {
        throw std::runtime_error("Could not open backbone torsion file.");
    }
    //bbTorsions = {};
    size_t res1, res2, res3;
    size_t plane_angle_index, bisect_angle_index, phi_index, psi_index;
    int energy, obs, tot_obs, val;
    while (triple_file >> res1 >> res2 >> res3 >> plane_angle_index >> bisect_angle_index >> phi_index >> psi_index >> energy >> obs >> tot_obs >> val) {
        // Store the value in the array
        if (res1 < 20 && res2 < 20 && res3 < 20 && plane_angle_index < 6 && bisect_angle_index < 6 && phi_index < 6 && psi_index < 6) {
            bbTorsions[res1][res2][res3][plane_angle_index][bisect_angle_index][phi_index][psi_index] = energy;
        } else {
            std::cerr << "Index out of bounds: " << res1 << " " << res2 << " " << res3 << " " << plane_angle_index << " " << bisect_angle_index << " " << phi_index << " " << psi_index << std::endl;
        }
    }
    triple_file.close();
    // Load sidechain torsion parameters
    std::ifstream sctorsion_file(sctorsion_filename);
    if (!sctorsion_file.is_open()) {
        throw std::runtime_error("Could not open sidechain torsion file.");
    }
    //scTorsions = {};
    size_t chi1, chi2, chi3, chi4;
    while (sctorsion_file >> res1 >> res2 >> res3 >> chi1 >> chi2 >> chi3 >> chi4 >> energy >> obs >> tot_obs >> val) {
        // Store the value in the array
        if (res1 < 20 && res2 < 20 && res3 < 20 && chi1 < 12 && chi2 < 12 && chi3 < 12 && chi4 < 12) {
            scTorsions[res1][res2][res3][chi1][chi2][chi3][chi4] = energy;
        } else {
            std::cerr << "Index out of bounds: " << res1 << " " << res2 << " " << res3 << " " << chi1 << " " << chi2 << " " << chi3 << " " << chi4 << std::endl;
        }
    }
    sctorsion_file.close();
}

void MCPUForceField::loadAromaticPotential(const std::string& filename) {
    // Load aromatic potential implementation
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open aromatic file.");
    }
    // Define variables to read the file
    int index, value, obs, tot_obs;
    // Read file
    while (file >> index >> value >> obs >> tot_obs) {
        aromaticPotential[index] = value;
    }
    file.close();
}

std::tuple<std::vector<int>, std::vector<int>> MCPUForceField::getSmogType(Topology& topology) {
    std::vector<int> remove_atom_ids, smogTypeVector;
    std::string sidechain_key, backbone_key;
    // Assign Smog Type
    for (int i = 0; i < topology.getNumAtoms(); ++i) {
        sidechain_key = topology.residues[topology.atoms[i].resID].resName + std::string("_") + topology.atoms[i].atomName;
        backbone_key = std::string("XXX_") + topology.atoms[i].atomName;
        if (smogTypeMap.find(sidechain_key) != smogTypeMap.end()) {
            int smogType = smogTypeMap[sidechain_key];
            smogTypeVector.push_back(smogType);
        } else if (smogTypeMap.find(backbone_key) != smogTypeMap.end()) {
            int smogType = smogTypeMap[backbone_key];
            smogTypeVector.push_back(smogType);
        } else {
            std::cerr << "Error: Smog Type not found for atom " << topology.atoms[i].atomID << std::endl;
            remove_atom_ids.push_back(topology.atoms[i].atomID);
        }
    }
    return std::make_tuple(smogTypeVector, remove_atom_ids);
}

std::tuple<Topology, Eigen::Matrix3Xd> MCPUForceField::removeAtomsByID(Topology& old_topology, Eigen::Matrix3Xd old_positions, std::vector<int> remove_atom_ids) {
    // Remove atoms by ID implementation
    Topology new_topology = old_topology;
    Eigen::Matrix3Xd new_positions;
    new_positions.resize(3, old_topology.getNumAtoms() - remove_atom_ids.size());
    std::vector<Atom> new_atoms;
    std::vector<Residue> new_residues;
    int new_atom_id = 0;
    std::unordered_map<int, int> old_to_new_atom_id;
    for (int i = 0; i < old_topology.getNumAtoms(); ++i) {
        int old_atom_id = old_topology.atoms[i].atomID;
        if (std::find(remove_atom_ids.begin(), remove_atom_ids.end(), old_topology.atoms[i].atomID) == remove_atom_ids.end()) {
            new_atoms.push_back(old_topology.atoms[i]);
            new_atoms.back().atomID = new_atom_id;
            new_positions.col(new_atom_id) = old_positions.col(i);
            old_to_new_atom_id[old_atom_id] = new_atom_id;
            new_atom_id++;
        }
    }
    new_topology.atoms = new_atoms;
    for (auto& residue : new_topology.residues) {
        for (auto& atom : residue.atoms) {
            atom.atomID = old_to_new_atom_id[atom.atomID];
        }
    }
    new_topology.updateTopology();
    return std::make_tuple(new_topology, new_positions);
}

Eigen::MatrixXd MCPUForceField::getMuPotential() const {
    return muPotential;
}

BBTorsionArray MCPUForceField::getBackBoneTorsions() const {
    return bbTorsions;
}

SCTorsionArray MCPUForceField::getSideChainTorsions() const {
    return scTorsions;
}

std::unordered_map<std::string, int> MCPUForceField::getSmogTypeMap() const {
    return smogTypeMap;
}

std::map<int, double> MCPUForceField::getHBondPotential() const {
    return hbondPotential;
}

HBondParamArray MCPUForceField::getSequenceDependentHBondPotential() const {
    return sequenceDependentHBondPotential;
}

std::map<int, double> MCPUForceField::getAromaticPotential() const {
    return aromaticPotential;
}

HBondConfig MCPUForceField::getHBondConfig() const {
    return hbondConfig;
}