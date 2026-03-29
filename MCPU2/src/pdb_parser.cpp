#include "pdb_parser.h"
#include "utils/string_utils.h"
#include "amino_acids.h"

#include <sstream>
#include <vector>
#include <iostream>
#include <fstream>
#include <set>
#include <tuple>

PDBFile::PDBFile(const std::string &filename) {
    parse_pdb_(filename);
    create_protein_structure_();
    create_positions_();
}

PDBFile::~PDBFile() {
    // Destructor implementation
}

void PDBFile::parse_pdb_(const std::string &filename) {
    // Read the PDB file and extract atom information
    // atom_number and residue_number are direct readouts from the PDB file
    // Indices used for the simulation will be 0-based and will be reassigned in create_topology_ function.
    std::ifstream pdbFile(filename);
    if (!pdbFile.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        exit(2);
    }
    // Read lines
    std::string line;
    while (std::getline(pdbFile, line)) {
        // Process only ATOM or HETATM records
        if (line.substr(0, 6) == "ATOM  " || line.substr(0, 6) == "HETATM") {
            // Extract atom information based on fixed-column positions in the PDB format
            int atom_number = std::stoi(line.substr(6, 5)); // Atom serial number (columns 7-11)
            std::string atom_name = trim(line.substr(12, 4));  // Atom name (columns 13-16)
            std::string element = trim(line.substr(76, 2));   // Element (columns 77-78)
            std::string res_name = trim(line.substr(17, 3));   // Residue name (columns 18-20)
            char chain_name = line[21];           // Chain ID (column 22)
            int residue_number = std::stoi(line.substr(22, 4)); // Residue sequence number (columns 23-26)
            Eigen::Vector3f coord;
            coord(0) = std::stod(line.substr(30, 8)); // X coordinate (columns 31-38)
            coord(1) = std::stod(line.substr(38, 8));      // Y coordinate (columns 39-46)
            coord(2) = std::stod(line.substr(46, 8));      // Z coordinate (columns 47-54)
            // Add atom information to the vector
            pdb_atoms_.push_back(PDBAtom(atom_number, atom_name, element, res_name, residue_number, chain_name, coord));
        }
    }
    // Generate a molecular tree with the parsed data
    // This will create a map of atoms 
    std::set<std::tuple<char, int>> residue_set;
    std::set<char> chain_set;
    for (const auto& pdb_atom : pdb_atoms_) {
        residue_set.insert({pdb_atom.chain_name, pdb_atom.residue_number});
        chain_set.insert(pdb_atom.chain_name);
    }
    // Map char chain name to chain ID
    int chain_id = 0;
    for (const auto& chain : chain_set) {
        chain_name_to_id_[chain] = chain_id++;
    }
    // Clear memory
    residue_set.clear();
    chain_set.clear();
    pdbFile.close();
}

void PDBFile::create_protein_structure_() {
    std::map<char, Chain*> chain_map;
    std::map<std::tuple<char, int>, Residue*> residue_map;

    for (const auto& pdb_atom : pdb_atoms_) {
        // Chain
        Chain* chain_ptr = nullptr;
        if (chain_map.find(pdb_atom.chain_name) == chain_map.end()) {
            auto chain = std::make_unique<Chain>(pdb_atom.chain_name, chain_name_to_id_[pdb_atom.chain_name]);
            chain_ptr = chain.get();
            protein.add_chain(std::move(chain));
            chain_map[pdb_atom.chain_name] = chain_ptr;
        } else {
            chain_ptr = chain_map[pdb_atom.chain_name];
        }

        // Residue
        Residue* residue_ptr = nullptr;
        auto residue_key = std::make_tuple(pdb_atom.chain_name, pdb_atom.residue_number);
        if (residue_map.find(residue_key) == residue_map.end()) {
            auto residue = std::make_unique<Residue>(pdb_atom.res_name, pdb_atom.residue_number, chain_name_to_id_[pdb_atom.chain_name]);
            residue_ptr = residue.get();
            chain_ptr->add_residue(std::move(residue));
            residue_map[residue_key] = residue_ptr;
        } else {
            residue_ptr = residue_map[residue_key];
        }

        // Atom
        auto atom = std::make_unique<Atom>(
            pdb_atom.atom_name,
            pdb_atom.element,
            pdb_atom.atom_number,
            pdb_atom.residue_number
        );
        residue_ptr->add_atom(std::move(atom));
    }
    protein.recompute_number_of_particles();
}

void PDBFile::create_positions_() {
    // Create positions from the parsed PDB data
    positions.resize(3, pdb_atoms_.size());
    for (size_t i = 0; i < pdb_atoms_.size(); ++i) {
        positions.col(i) = pdb_atoms_[i].coord;
    }
}