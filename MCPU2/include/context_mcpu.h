#ifndef CONTEXT_MCPU_H
#define CONTEXT_MCPU_H

#include <Eigen/Dense>
#include <vector>

#include "cell_list.h"
#include "system_mcpu.h"

struct Atom;
struct Residue;

struct Energy {
    float mu_energy;
    float bb_torsion_energy;
    float sc_torsion_energy;
    float aromatic_energy;
    float hbond_energy;
    float total_energy;
};

struct HBondState {
    // Illustration can be found at https://doi.org/10.1016/j.str.2006.11.010
    bool                  hbond_state;
    size_t                donor_resType;
    size_t                acceptor_resType;
    int                   secondary_structure_label;
    std::array<double, 3> bisect_angles;
    std::array<double, 3> plane_angles;
};

class ContextMCPU {
   public:
    ContextMCPU(SystemMCPU &system, const Eigen::Matrix3Xf positions);
    ~ContextMCPU();

    // Context of the system
    std::vector<std::vector<Eigen::Vector3f>> backbone_plane_vectors_per_chain;   // N, CA, O
    std::vector<std::vector<Eigen::Vector3f>> backbone_bisect_vectors_per_chain;  // N, CA, O
    std::vector<std::vector<Eigen::Vector3f>>
        plane_vectors_per_chain;  // N, CA, C for hydrogen bond
    std::vector<std::vector<Eigen::Vector3f>>
        bisect_vectors_per_chain;  // N, CA, C for hydrogen bond

    // Pairwise interaction features
    Eigen::MatrixXf pairwise_contact_matrix;
    Eigen::MatrixXf ca_ca_distance2_matrix;
    Eigen::MatrixXi hydrogen_bond_feature_matrix;
    // Backbone torsion features
    std::vector<std::vector<float>> backbone_phi_angles_per_chain;
    std::vector<std::vector<float>> backbone_psi_angles_per_chain;
    std::vector<std::vector<int>>   backbone_phi_feats_per_chain;
    std::vector<std::vector<int>>   backbone_psi_feats_per_chain;
    std::vector<std::vector<int>>   backbone_plane_feats_per_chain;
    std::vector<std::vector<int>>   backbone_bisect_feats_per_chain;
    // Sidechain torsion features
    std::vector<std::vector<std::array<int, 4>>> sidechain_torsion_feats_per_chain;
    // Aromatic interaction features
    std::vector<Eigen::Vector3f> aromatic_centers;
    std::vector<int>             aromatic_interaction_feats;
    // Hydrogen bond features
    std::vector<HBondState> hydrogen_bond_feats;

    Energy energy;

    void initialize();
    void update();
    void reorder_positions();
    void compute_pairwise_contacts(int atom_id);
    void compute_backbone_vectors(int chain_id, int residue_id);
    void compute_vectors(int chain_id, int residue_id);
    void compute_backbone_angles(int chain_id, int residue_id);
    void compute_sidechain_angles(int chain_id, int residue_id);
    void compute_hydrogen_bonds(int chain_id, int residue_id);
    void compute_aromatic_centers();
    void compute_aromatic_angles();

    std::vector<int> filter_out_end_of_chain_residues(std::vector<int> residue_ids, int n_residues);
    std::vector<int> filter_out_c_terminal_residue(std::vector<int> residue_ids, int n_residues);
    std::vector<int> filter_out_n_terminal_residue(std::vector<int> residue_ids, int n_residues);

    void compute_energy(int chain_id, std::vector<int> moved_residue_ids);

   private:
    SystemMCPU      &system_;
    Eigen::Matrix3Xf positions_;
    CellList         cell_list_;
    int              n_chains_;
    int              n_residues_;
    int              n_atoms_;
};

#endif  // CONTEXT_MCPU_H
