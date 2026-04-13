#include "db/amino_acid_db.h"
#include <fstream>
#include <iostream>

namespace mcpu::db {

std::vector<AminoAcidDB> LoadTopology(const std::string& filepath) {
    std::ifstream file(filepath);
    
    // Safety check for cluster paths
    if (!file.is_open()) {
        std::cerr << "ERROR: Could not find topology file at: " << filepath << std::endl;
        return {}; 
    }

    try {
        nlohmann::json j;
        file >> j;

        // If your JSON has a root key "amino_acids", use j.at("amino_acids").get<...>()
        // Otherwise, if the JSON is just a list [...], use this:
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

}