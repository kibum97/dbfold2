#ifndef LOAD_PDB_H
#define LOAD_PDB_H

#include <Eigen/Dense>
#include <map>
#include <string>

#include "protein_structure.h"

struct PDBAtom {
    int             atom_number;
    std::string     atom_name, element, res_name;
    int             residue_number;
    char            chain_name;
    Eigen::Vector3f coord;

    PDBAtom(int atom_number_, const std::string &atom_name_, const std::string &element_,
            const std::string &res_name_, int residue_number_, const char &chain_name_,
            const Eigen::Vector3f &coord_)
        : atom_number(atom_number_),
          atom_name(atom_name_),
          element(element_),
          res_name(res_name_),
          residue_number(residue_number_),
          chain_name(chain_name_),
          coord(coord_) {}
};  // Structure to hold atom information from PDB file

class PDBFile {
   public:
    PDBFile(const std::string &filename);
    ~PDBFile();

    ProteinStructure protein;    // Topology object to hold the parsed topology
    Eigen::Matrix3Xf positions;  // Matrix to hold the positions of atoms

   private:
    std::map<char, int> chain_name_to_id_;  // Map to convert chain name to ID

    void                 parse_pdb_(const std::string &filename);
    void                 create_protein_structure_();
    void                 create_positions_();
    std::vector<PDBAtom> pdb_atoms_;  // Vector to hold atom information
};

#endif  // LOAD_PDB_H