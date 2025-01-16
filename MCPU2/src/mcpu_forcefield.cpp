#include "mcpu_forcefield.h"
#include "read_config.h"
#include <fstream>
#include <iostream>

//TODO: Need to read hbond_file from config file - Need to add this var

int USE_SMOG_TYPE = 1;
extern std::string potential_file;
extern std::string atom_type_file;
extern std::string hbond_file;
extern std::string hydrogen_bonding_data;

MCPUForceField::MCPUForceField() {
    // Constructor implementation
}

MCPUForceField::~MCPUForceField() {
    // Destructor implementation
}

void MCPUForceField::initialize() {
    // Initialize force field implementation
    loadMuPotential(potential_file);
    loadSmogTypeMap(atom_type_file);
    loadHbondPotential(hbond_file);
}

void MCPUForceField::loadMuPotential(const std::string& filename) {
    // Load mu potential implementation
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file.");
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
        throw std::runtime_error("Could not open file.");
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
        throw std::runtime_error("Could not open file.");
    }

    // Find the matrix dimensions dynamically
    std::map <int, double> hbondPotential;
    int ss;
    int phi_donor, psi_donor, phi_acceptor, psi_acceptor;
    int bisect_angle, plane_angle;
    double value;
    int index;

    // Store the value in the map
    while (file >> ss >> phi_donor >> psi_donor >> phi_acceptor >> psi_acceptor >> bisect_angle >> plane_angle >> value) {    
        index = computeHBondindex(ss, phi_donor, psi_donor, phi_acceptor, psi_acceptor, bisect_angle, plane_angle);
        hbondPotential[index] = value;
    }
    file.close();    
}

int MCPUForceField::computeHBondindex(int ss, int phi_donor, int psi_donor, int phi_acceptor, int psi_acceptor, int bisect_angle, int plane_angle) {
    // Compute hydrogen bond index implementation
    return plane_angle 
       + bisect_angle * 9 
       + psi_acceptor * 9 * 9 
       + phi_acceptor * 9 * 9 * 9 
       + psi_donor * 9 * 9 * 9 * 9 
       + phi_donor * 9 * 9 * 9 * 9 * 9 
       + ss * 9 * 9 * 9 * 9 * 9 * 9;
}

Eigen::MatrixXd MCPUForceField::getMuPotential() const {
    return muPotential;
}

std::unordered_map<std::string, int> MCPUForceField::getSmogTypeMap() const {
    return smogTypeMap;
}

std::map<int, double> MCPUForceField::getHBondPotential() const {
    return hbondPotential;
}

HBondConfig MCPUForceField::getHBondConfig() const {
    return hbondConfig;
}