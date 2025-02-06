#ifndef ROTAMER_H
#define ROTAMER_H

#include <string>
#include <array>
#include <vector>
#include <iostream>
#include <unordered_map>
#include <fstream>
#include <sstream>

struct RotamerData {
    double probability;
    std::array<double, 4> chiMeans;
    std::array<double, 4> chiStdDevs;
};

//const std::unordered_map<std::string, std::vector<RotamerData>> readRotamerData(const std::string& rotamerFile);

inline const std::unordered_map<std::string, std::vector<RotamerData>> readRotamerData(const std::string& rotamerFile) {
    // Read rotamer data from file
    std::unordered_map<std::string, std::vector<RotamerData>> rotamerMap;
    std::ifstream file(rotamerFile);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << rotamerFile << " for reading." << std::endl;
        return rotamerMap;
    }
    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string resName;
        int rotLib, rotType, angleSet, angleIdx, nObs, nOcc; // This will be ignored
        double prob, chi1, stddev1, chi2, stddev2, chi3, stddev3, chi4, stddev4;
        if (iss >> resName >> rotLib >> rotType >> angleSet >> angleIdx >> nObs >> nOcc >> prob >> chi1 >> stddev1 >> chi2 >> stddev2 >> chi3 >> stddev3 >> chi4 >> stddev4) {
            RotamerData rotamerData = {prob, {chi1, chi2, chi3, chi4}, {stddev1, stddev2, stddev3, stddev4}};
            rotamerMap[resName].push_back(rotamerData);
        }
    }
    return rotamerMap;
}

#endif // ROTAMER_H