#include "protein_structure.h"
#include <iostream>

void ProteinStructure::add_chain(std::unique_ptr<Chain> chain) {
    chains.emplace_back(std::move(chain));
}

void ProteinStructure::add_residue_to_chain(char chain_id, std::unique_ptr<Residue> residue) {
    for (const auto& chain : chains) {
        if (chain->id == chain_id) {
            chain->add_residue(std::move(residue));
            return;
        }
    }
}

void ProteinStructure::add_atom_to_residue(char chain_id, int residue_id, std::unique_ptr<Atom> atom) {
    for (const auto& chain : chains) {
        if (chain->id == chain_id) {
            for (const auto& residue : chain->residues) {
                if (residue->id == residue_id) {
                    residue->add_atom(std::move(atom));
                    return;
                }
            }
        }
    }
}

void ProteinStructure::validate_structure() const {
    // Check if the structure is correctly formatted
    // Each chain should have full residues without missing residues
    // Residue IDs will be 0-based during this check
    int start_residue_id = 0;
    int start_atom_id = chains[0]->residues[0]->atoms[0]->id;
    for (auto& chain: chains) {
        start_residue_id = chain->residues[0]->id;
        for (size_t i = 0; i < chain->residues.size(); ++i) {
            if (chain->residues[i]->id != start_residue_id + i) {
                std::cerr << "Error: Missing residue in the chain (Residue ID mismatch): chain " << chain->id << std::endl;
                exit(EXIT_FAILURE);
            }
            // Check if the residue ID is unique
            for (size_t j = i + 1; j < chain->residues.size(); ++j) {
                if (chain->residues[i]->id == chain->residues[j]->id) {
                    std::cerr << "Error: Duplicate residue ID in the chain: chain " << chain->id << std::endl;
                    exit(EXIT_FAILURE);
                }
            }
        }
        // Each residue should have full atoms without missing atoms
        // Atom IDs should be unique (globally)
        // Atom IDs will be 0-based during this check
        size_t global_atom_id = start_atom_id;
        for (auto& residue: chain->residues) {
            for (size_t k = 0; k < residue->atoms.size(); ++k) {
                if (residue->atoms[k]->id != global_atom_id) {
                    std::cerr << "Error: Missing atom in the residue (Atom ID mismatch) chain " << chain->id << " residue " << residue->id << std::endl;
                    exit(EXIT_FAILURE);
                }
                // Check if the atom ID is unique
                for (size_t l = k + 1; l < residue->atoms.size(); ++l) {
                    if (residue->atoms[k]->id == residue->atoms[l]->id) {
                        std::cerr << "Error: Duplicate atom ID in the residue: chain " << chain->id << " residue " << residue->id << std::endl;
                        exit(EXIT_FAILURE);
                    }
                }
                global_atom_id++;
            }
        }
        // Check if the chain ends with a OXT atom
        if (!chain->residues.empty()) {
            const auto& last_residue = chain->residues.back();
            bool ends_with_oxt = false;
            for (const auto& atom : last_residue->atoms) {
                if (atom->name == "OXT") {
                    ends_with_oxt = true;
                    break;
                }
            }
            if (!ends_with_oxt) {
                std::cerr << "Error: Last residue in the chain does not have OXT atom: chain " << chain->id << std::endl;
                exit(EXIT_FAILURE);
            }
        }
    }
    // If the structure is valid, make residue IDs and atom IDs 0-based
    for (auto& chain: chains) {
        // Residue IDs are 0-based per chain
        // Atom IDs are 0-based globally
        start_residue_id = chain->residues[0]->id;
        for (auto& residue: chain->residues) {
            residue->id -= start_residue_id;
            for (auto& atom: residue->atoms) {
                atom->id -= start_atom_id;
            }
        }
    }
}

void ProteinStructure::recompute_number_of_particles() {
    num_atoms = 0;
    num_residues = 0;
    num_chains = 0;
    for (const auto& chain : chains) {
        int num_atoms_of_chain = 0;
        num_residues_per_chain.push_back(chain->residues.size());
        num_residues += chain->residues.size();
        for (const auto& residue : chain->residues) {
            num_atoms_of_chain += residue->atoms.size();
            num_atoms += residue->atoms.size();
        }
        num_atoms_per_chain.push_back(num_atoms_of_chain);
        num_chains++;
    }
}

void ProteinStructure::print_structure() const {
    for (const auto& chain : chains) {
        std::cout << "Chain " << chain->id << ":\n";
        for (const auto& residue : chain->residues) {
            std::cout << "  Residue " << residue->name << " " << residue->id << ":\n";
            for (const auto& atom : residue->atoms) {
                std::cout << "    Atom " << atom->name << " (" << atom->element << ") ID: " << atom->id << "\n";
            }
        }
    }
    std::cout << "Number of chains: " << num_chains << "\n";
    std::cout << "Number of residues: " << num_residues << "\n";
    std::cout << "Number of atoms: " << num_atoms << "\n";
    std::cout << "Number of residues per chain: ";
    for (const auto& num_residues : num_residues_per_chain) {
        std::cout << num_residues << " ";
    }
    std::cout << "\n";
    std::cout << "Number of atoms per chain: ";
    for (const auto& num_atoms : num_atoms_per_chain) {
        std::cout << num_atoms << " ";
    }
    std::cout << "\n";
}