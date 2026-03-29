#ifndef SYSTEM_MCPU_H
#define SYSTEM_MCPU_H

#include <Eigen/Dense>
#include <map>

#include "forces/ForcefieldMCPU.h"

struct Atom;
struct Residue;
class ProteinStructure;

extern struct HBondConfig {
    const std::vector<std::string> donors          = {"N", "CA", "C", "C", "CA", "N", "CA"};
    const std::vector<std::string> acceptors       = {"O", "C", "CA", "N", "N", "CA", "C", "CA"};
    const std::vector<int>         donor_offset    = {0, 0, -1, 0, -1, -1, 1};
    const std::vector<int>         acceptor_offset = {0, 0, 0, 1, 0, 1, 1, -1};
    const float                    min_D           = 0.0;
    const float                    max_D           = 5.0;
    const float                    D_int           = 0.1;
} HBOND_CONFIG;

struct mcpu_atom {
    int   chain_id;      // Chain ID
    int   residue_id;    // global residue id; residues in later chain will have shifted residue ids
    int   col_index;     // column index in the positions matrix
    Atom *atom;          // Pointer to Atom object
    Residue *residue;    // Pointer to Residue object
    int      smog_type;  // Smog type
    float    radius;     // Atom radius

    // Default constructor
    mcpu_atom()
        : chain_id(-1),
          residue_id(-1),
          col_index(-1),
          atom(nullptr),
          residue(nullptr),
          smog_type(-1),
          radius(0.0f) {}
    // Constructor to initialize the mcpu_atom structure
    mcpu_atom(int chain_id_, int residue_id_, int col_index_, Atom *atom_, Residue *residue_,
              int smog_type_, float radius_)
        : chain_id(chain_id_),
          residue_id(residue_id_),
          col_index(col_index_),
          atom(atom_),
          residue(residue_),
          smog_type(smog_type_),
          radius(radius_) {}
};

struct mcpu_hydrogen_bond {
    int                chain_id;                   // Chain ID to read values from context
    int                residue_id;                 // Residue ID to read values from context
    int                acceptor_o_id;              // ID of the acceptor oxygen atom
    std::array<int, 3> donor_n_triangle_atoms;     // ID of atoms forming a triangle with the donor
                                                   // nitrogen as the apex
    std::array<int, 3> acceptor_o_triangle_atoms;  // ID of atoms forming a triangle with the
                                                   // acceptor oxygen as the apex
    // Other variables needed can be accessed using residue_id
};

struct mcpu_dihedral_angles {
    std::array<int, 4> dihedral_atom_ids;
};

struct mcpu_aromatic_atoms {
    std::array<int, 3> aromatic_atom_ids;
};

class SystemMCPU {
   public:
    SystemMCPU(ForcefieldMCPU &forcefield, ProteinStructure &protein_structure);
    ~SystemMCPU();

    // Member functions
    // Init functions
    void pdb_atoms_to_mcpu_atoms();
    void init_backbone_atoms();
    void init_sidechain_dihedrals();
    void init_aromatic_interactions();
    void init_hydrogen_bonds();
    void init_system_forces();
    void init_mu_potential_coefficients();
    void init_torsion_parameters();
    void init_aromatic_parameters();
    void init_hydrogen_bond_parameters();

    // Util functions
    int get_global_residue_index(int chain_id, int residue_id);
    int get_global_residue_id(int chain_id, int residue_id);

    // Member variables
    // Vectors to store index for positions matrix
    // These are global for all chains and pair of start and end indices
    std::vector<int> chain_indices;    // List of chain indices
    std::vector<int> residue_indices;  // List of residue indices
    // Vectors to store atom information
    std::vector<mcpu_atom> mcpu_atoms;  // List of atoms that are read from the PDB file
    // Vectors to store atom ids for each feature
    std::vector<std::vector<int>>
        backbone_atoms;  // List of backbone N-CA-C atoms per chain; plane/bisect vectors will be
                         // calculated based on these atoms
    std::vector<std::vector<int>>
        sidechain_atom_indices;  // List of first sidechain atom indices of residues per chain
    std::vector<std::vector<int>>
        number_of_sidechain_dihedrals_per_chain;  // Number of sidechain dihedrals per chain
    std::vector<mcpu_aromatic_atoms> aromatic_interactions;  // List of aromatic atoms
    std::vector<mcpu_hydrogen_bond>  hydrogen_bonds;         // List of hydrogen bond atoms

    // Number of particles per chain
    std::vector<int> n_residues_per_chain;  // Number of residues per chain
    std::vector<int> n_atoms_per_chain;     // Number of atoms per chain
    // PDB write order - have to think about whether this is needed
    std::map<int, int>
        pdb_atom_id_to_system_atom_index;  // Map to convert PDB atom ID to system atom index
    // Array to store contact cutoff2 matrix
    ContactCutoff2Matrix contact_cutoff2;  // Contact cutoff matrix for pairwise interactions
    // Force field reference
    // These are formatted to be compatible with the ContextMCPU class
    Eigen::MatrixXf mu_potential_coefficients;  // Mu potential coefficients
    std::vector<std::vector<BackboneTorsionPerTriplet>>  backbone_torsion_parameters_per_chain;
    std::vector<std::vector<SidechainTorsionPerTriplet>> sidechain_torsion_parameters_per_chain;
    std::vector<std::vector<int>>                        residue_types_per_chain;
    std::map<int, float>                                 hydrogen_bond_parameters;
    HBondParams                                          sequence_based_hydorgen_bond_parameters;
    std::map<int, float>                                 aromatic_parameters;

   private:
    ForcefieldMCPU   &forcefield_;         // Reference to the force field
    ProteinStructure &protein_structure_;  // Reference to the protein structure
    int               n_chains_;           // Number of chains
    int               n_residues_;         // Number of residues
    int               n_atoms_;            // Number of atoms
};

#endif  // SYSTEM_MCPU_H