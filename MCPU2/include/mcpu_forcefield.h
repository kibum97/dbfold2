#ifndef FORCEFIELD_H
#define FORCEFIELD_H

#include <Eigen/Dense>
#include <string>
#include <map>
#include "topology.h"

extern int USE_SMOG_TYPE;

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
    std::unordered_map<std::string, int> getSmogTypeMap() const;
    std::tuple<std::vector<int>, std::vector<int>> getSmogType(Topology& topology);

    // Updating topology and positions
    std::tuple<Topology, Eigen::Matrix3Xd> removeAtomsByID(Topology& old_topology, Eigen::Matrix3Xd old_positions, std::vector<int> remove_atom_ids);

    std::map<int, double> getHBondPotential() const;
    HBondConfig getHBondConfig() const;
    int computeHBondindex(int ss, int phi_donor, int psi_donor, int phi_acceptor, int psi_acceptor, int bisect_angle, int plane_angle);

private:
    //Eigen::VectorXd bbTorsions;
    //Eigen::VectorXd scTorsions;
    //Eigen::VectorXi aromatics; // 0 if residue is not aromatic, 1 if residue is aromatic
    //Eigen::VectorXd hbonds;
    Eigen::MatrixXd muPotential;
    std::unordered_map<std::string, int> smogTypeMap;
    struct HBondConfig hbondConfig;
    std::map<int, double> hbondPotential;

    // Contact Energy - Mu Potential
    void loadMuPotential(const std::string& filename);
    void loadSmogTypeMap(const std::string& filename);
    // Hydrogen Bonding Energy
    void loadHbondPotential(const std::string& filename);

    //void loadAromatics(const std::string& filename);
};

#endif // FORCEFIELD_H