#include "mcpu_aromatic.h"
#include "topology.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cstdlib>

std::unordered_map<int, int> InitAromatic(const std::string& filename) {
    std::unordered_map<int, int> aromatic_energy;
    std::ifstream file(filename);
    std::string line;

    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        std::exit(EXIT_FAILURE);
    }

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        int plane_angle_bin, value;
        if (iss >> plane_angle_bin >> value) {
            aromatic_energy[plane_angle_bin] = value;
        }
    }

    file.close();
    return aromatic_energy;
}

std::vector<int> AromaticsFromTopology(const std::vector<Atom>& atoms) {
    std::vector<int> aromatics(atoms.size());
    for (size_t i = 0; i < atoms.size(); ++i) {
        if (atoms[i].resName == "PHE" || atoms[i].resName == "TYR" || atoms[i].resName == "TRP") {
            aromatics[i] = 1;
        } else {
            aromatics[i] = 0;
        }
    }
    return aromatics;
}

 Eigen::VectorXi computeStepFunction(const Eigen::VectorXd& values, double stepSize, double offset = 0.0) {
    return ((values.array() - offset) / stepSize).floor().cast<int>();
}