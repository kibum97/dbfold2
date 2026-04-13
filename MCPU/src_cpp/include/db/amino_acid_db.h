#pragma once
#include <nlohmann/json.hpp>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>

namespace mcpu::db {

struct Torsion {
    int id;
    std::vector<std::string> atoms;
    std::vector<std::string> affected_atoms;
};

struct AminoAcidDB {
    std::string name;
    std::string symbol;
    int heavy_atoms;
    int ntorsions;
    std::vector<Torsion> torsions;
};

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Torsion, id, atoms, affected_atoms)
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(AminoAcidDB, name, symbol, heavy_atoms, ntorsions, torsions)

inline std::vector<AminoAcidDB> LoadTopology(const std::string& filepath) {
    std::ifstream file(filepath);
    
    // Safety check for cluster paths
    if (!file.is_open()) {
        std::cerr << "ERROR: Could not find topology file at: " << filepath << std::endl;
        return {}; 
    }

    try {
        nlohmann::json j;
        file >> j;

        if (j.is_array()) {
            return j.get<std::vector<AminoAcidDB>>();
        } else if (j.contains("amino_acids")) {
            return j.at("amino_acids").get<std::vector<AminoAcidDB>>();
        }
        
        return {};
    } catch (const nlohmann::json::exception& e) {
        std::cerr << "JSON Parse Error in " << filepath << ": " << e.what() << std::endl;
        return {};
    }
}

} // namespace mcpu::db