#include "forces/ForcefieldMCPU.h"
//#include "protein_structure.h"
#include "utils/mcpu_hbond_utils.h"
#include <fstream>
#include <iostream>

extern std::string potential_file;
extern std::string atom_type_file;
extern std::string hbond_file; 
extern std::string seq_hbond_file;
extern std::string hydrogen_bonding_data;
extern std::string all_triplet_file;
extern std::string all_sctorsion_file;

ForcefieldMCPU::ForcefieldMCPU() {
    // Initialize the torsion arrays
    std::fill(&bbTorsions[0][0][0][0][0][0][0], &bbTorsions[0][0][0][0][0][0][0] + sizeof(bbTorsions) / sizeof(bbTorsions[0][0][0][0][0][0][0]), 1000);
    std::fill(&scTorsions[0][0][0][0][0][0][0], &scTorsions[0][0][0][0][0][0][0] + sizeof(scTorsions) / sizeof(scTorsions[0][0][0][0][0][0][0]), 1000);    
    // Initialize force field implementation
    loadMuPotential(potential_file);
    loadSmogTypeMap(atom_type_file);
    initialize_contact_cutoff2();
    loadHbondPotential(hbond_file);
    loadSequenceDependentHBondPotential(seq_hbond_file);
    loadTorsionalPotential(all_triplet_file, all_sctorsion_file);
}

ForcefieldMCPU::~ForcefieldMCPU() {
    // Destructor implementation
}

void ForcefieldMCPU::loadMuPotential(const std::string& filename) {
    // Load mu potential implementation
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open mu potential file.");
    }

    // Find the matrix dimensions dynamically
    std::vector<std::tuple<int, int, float>> data;
    int maxRow = 0, maxCol = 0;
    int row, col;
    float value;

    while (file >> row >> col >> value) {
        data.emplace_back(row, col, value);
        if (row > maxRow) maxRow = row;
        if (col > maxCol) maxCol = col;
    }
    file.close();

    // Create Eigen matrix with the appropriate size
    Eigen::MatrixXf muParamMatrix = Eigen::MatrixXf::Zero(maxRow + 1, maxCol + 1);

    // Populate the matrix
    for (const auto& [rowIdx, colIdx, val] : data) {
        muParamMatrix(rowIdx, colIdx) = val;
        muParamMatrix(colIdx, rowIdx) = val;
    }
    muPotential = muParamMatrix;
}

void ForcefieldMCPU::loadSmogTypeMap(const std::string& filename) {
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
        float radius;

        // Parse the three columns
        if (ss >> atomName >> residueName >> smogType >> radius) {
            // Concatenate the first wo columns to fdonor_orm the key
            std::string key = residueName + "_" + atomName;
            // Insert into the map
            std::pair<int, float> smogTypeData = std::make_pair(smogType, radius);
            smogTypeMap[key] = smogTypeData;
        } else {
            std::cerr << "Error: Malformed line: " << line << std::endl;
        }
    }

    // Close the file
    inputFile.close();
}

void ForcefieldMCPU::initialize_contact_cutoff2() {
    // Initialize the contact cutoff matrix
    contact_cutoff2.setZero();
    for (const auto& [key1, smog_type1] : smogTypeMap) {
        int smog_type_index1 = smog_type1.first;
        float smog_radius1 = smog_type1.second;
        for (const auto& [key2, smog_type2] : smogTypeMap) {
            int smog_type_index2 = smog_type2.first;
            float smog_radius2 = smog_type2.second;
            // Set the cutoff value in the matrix
            if (smog_type_index1 > NUM_SMOG_TYPES || smog_type_index2 > NUM_SMOG_TYPES) {
                std::cerr << "Error: SMOG type index out of bounds: " << smog_type_index1 << ", " << smog_type_index2 << std::endl;
                exit(EXIT_FAILURE);
            }
            for (int i = 0; i < NUM_SMOG_TYPES; ++i) {
                contact_cutoff2(smog_type_index1, smog_type_index2) = 0.75 * 1.8 * 0.75 * 1.8 * (smog_radius1 + smog_radius2) * (smog_radius1 + smog_radius2);
                contact_cutoff2(smog_type_index2, smog_type_index1) = contact_cutoff2(smog_type_index1, smog_type_index2); // Symmetric
            }
            // Ignore backbone-backbone interactions by setting cutoff to 0.
            if (smog_type_index1 > 78 && smog_type_index2 > 78) {
                contact_cutoff2(smog_type_index1, smog_type_index2) = 0.0;
                contact_cutoff2(smog_type_index2, smog_type_index1) = 0.0;
            }
        }
    }
}

void ForcefieldMCPU::loadHbondPotential(const std::string& filename) {
    // Load hydrogen bond potential implementation
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open Hbond file.");
    }
    // Find the matrix dimensions dynamically
    int ss;
    int plane_angle_ca, bisect_angle_ca, plane_angle_ca_off, bisect_angle_ca_off;
    int bisect_angle_strand, plane_angle_strand, obs, tot_obs;
    float value;
    int index;

    // Store the value in the map
    while (file >> ss >> plane_angle_ca >> bisect_angle_ca >> plane_angle_ca_off >> bisect_angle_ca_off >> plane_angle_strand >> bisect_angle_strand >> value >> obs >> tot_obs) {    
        std::array <int, 3> bisect_indices = {bisect_angle_ca, bisect_angle_ca_off, bisect_angle_strand};
        std::array <int, 3> plane_indices = {plane_angle_ca, plane_angle_ca_off, plane_angle_strand};
        index = compute_hbond_index(ss, bisect_indices, plane_indices);
        hbondPotential[index] = value;
    }
    file.close();
}

void ForcefieldMCPU::loadSequenceDependentHBondPotential(const std::string& filename) {
    // Load sequence dependent hydrogen bond potential implementation
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open sequence dependent Hbond file.");
    }
    // Read and store values
    int secondary_structure_label, donor_resType, acceptor_resType;
    float value;
    int obs, tot_obs;
    while (file >> secondary_structure_label >> donor_resType >> acceptor_resType >> value >> obs >> tot_obs) {
        sequenceDependentHBondPotential[secondary_structure_label][donor_resType][acceptor_resType] = (-1.) * value;
    }
    file.close();
}

void ForcefieldMCPU::loadTorsionalPotential(const std::string& triple_filename, const std::string& sctorsion_filename) {
    // Load backbone torsion parameters
    std::ifstream triple_file(triple_filename);
    if (!triple_file.is_open()) {
        throw std::runtime_error("Could not open backbone torsion file.");
    }
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

void ForcefieldMCPU::loadAromaticPotential(const std::string& filename) {
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