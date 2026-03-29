#ifndef PROTEIN_STRUCTURE_H
#define PROTEIN_STRUCTURE_H

#include <memory>
#include <string>
#include <vector>

struct Atom {
    std::string name;
    std::string element;
    int         id;
    int         residue_id;

    Atom(const std::string &name_, const std::string &element_, int id_, int residue_id_)
        : name(name_), element(element_), id(id_), residue_id(residue_id_) {}
};

struct Residue {
    std::string                        name;
    int                                id;
    int                                chain_id;
    std::vector<std::unique_ptr<Atom>> atoms;

    Residue(const std::string &name_, int id_, int chain_id_)
        : name(name_), id(id_), chain_id(chain_id_) {}

    void add_atom(std::unique_ptr<Atom> atom) { atoms.emplace_back(std::move(atom)); }
};

struct Chain {
    char                                  name;
    int                                   id;
    std::vector<std::unique_ptr<Residue>> residues;

    Chain(char name_, int id_) : name(name_), id(id_) {}

    void add_residue(std::unique_ptr<Residue> residue) {
        residues.emplace_back(std::move(residue));
    }
};

class ProteinStructure {
   public:
    std::vector<std::unique_ptr<Chain>> chains;

    // Member variables
    int              num_atoms;
    int              num_residues;
    int              num_chains;
    std::vector<int> num_atoms_per_chain;     // Vector to hold the number of atoms per chain
    std::vector<int> num_residues_per_chain;  // Vector to hold the number of residues per chain

    // Utility functions
    void add_chain(std::unique_ptr<Chain> chain);
    void add_residue_to_chain(char chain_id, std::unique_ptr<Residue> residue);
    void add_atom_to_residue(char chain_id, int residue_id, std::unique_ptr<Atom> atom);
    void recompute_number_of_particles();
    void validate_structure() const;
    void print_structure() const;
};

#endif  // PROTEIN_STRUCTURE_H