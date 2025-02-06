#ifndef FORCEFIELD_H
#define FORCEFIELD_H

#include <Eigen/Dense>
#include <string>
#include <map>
#include "topology.h"

extern int USE_SMOG_TYPE;
using BBTorsionParamArray = std::array<std::array<std::array<std::array<int, 6>, 6>, 6>, 6>;
using BBTorsionArray = std::array<std::array<std::array<BBTorsionParamArray, 20>, 20>, 20>;
using SCTorsionParamArray = std::array<std::array<std::array<std::array<int, 12>, 12>, 12>, 12>;
using SCTorsionArray = std::array<std::array<std::array<SCTorsionParamArray, 20>, 20>, 20>;
using HBondParamArray = std::array<std::array<std::array<int, 3>, 20>, 20>;

struct HBondConfig {
    const std::vector<std::string> donors = {"N", "CA", "C", "C", "CA", "N", "CA"};
    const std::vector<int> donor_offset = {0, 0, -1, 0, -1, -1, 1};
    const std::vector<std::string> acceptors = {"O", "C", "CA", "N", "N", "CA", "C", "CA"};
    const std::vector<int> acceptor_offset = {0, 0, 0, 1, 0, 1, 1, -1};
    const float min_D = 0.0;
    const float max_D = 5.0;
    const float D_int = 0.1;
};

class MCPUForceField {
public:
    MCPUForceField();
    ~MCPUForceField();

    Eigen::MatrixXd getMuPotential() const;
    BBTorsionArray getBackBoneTorsions() const;
    SCTorsionArray getSideChainTorsions() const;
    std::unordered_map<std::string, int> getSmogTypeMap() const;
    std::tuple<std::vector<int>, std::vector<int>> getSmogType(Topology& topology);

    // Updating topology and positions
    std::tuple<Topology, Eigen::Matrix3Xd> removeAtomsByID(Topology& old_topology, Eigen::Matrix3Xd old_positions, std::vector<int> remove_atom_ids);    

    std::map<int, double> getHBondPotential() const;
    HBondParamArray getSequenceDependentHBondPotential() const;
    std::map<int, double> getAromaticPotential() const;
    HBondConfig getHBondConfig() const;

private:
    //Eigen::VectorXi aromatics; // 0 if residue is not aromatic, 1 if residue is aromatic
    //Eigen::VectorXd hbonds;
    Eigen::MatrixXd muPotential;
    BBTorsionArray bbTorsions;
    SCTorsionArray scTorsions;
    std::unordered_map<std::string, int> smogTypeMap;
    struct HBondConfig hbondConfig;
    std::map<int, double> hbondPotential;
    HBondParamArray sequenceDependentHBondPotential;
    std::map<int, double> aromaticPotential;

    // Contact Energy - Mu Potential
    void loadMuPotential(const std::string& filename);
    void loadSmogTypeMap(const std::string& filename);
    // Hydrogen Bonding Energy
    void loadHbondPotential(const std::string& filename);
    void loadSequenceDependentHBondPotential(const std::string& filename);
    // Torsion Energy
    void loadTorsionalPotential(const std::string& triplet_file, const std::string& sctorsion_file);
    // Aromatic Energy
    void loadAromaticPotential(const std::string& filename);

    //void loadAromatics(const std::string& filename);
};

#endif // FORCEFIELD_H