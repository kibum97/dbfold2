#ifndef CONTEXT_H
#define CONTEXT_H

#include <vector>
#include <cmath>
#include <chrono>
#include <Eigen/Dense>
#include <unordered_set>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include "topology.h"
#include "cell_list.h"
#include "system.h"

using BackBoneAtoms = std::array<size_t, 4>;
using DihedralAtoms = std::array<size_t, 4>;

struct HBondState {
    // Illustration can be found at https://doi.org/10.1016/j.str.2006.11.010
    bool hbond_state;
    size_t donor_resType;
    size_t acceptor_resType;
    int secondary_structure_label;
    std::array<double, 3> bisect_angles;
    std::array<double, 3> plane_angles;
};

struct Energy {
    double mu_energy;
    double bb_torsion_energy;
    double sc_torsion_energy;
    double aromatic_energy;
    double hbond_energy;
    double total_energy;
};

class Context {
public:
    Context(Topology& topology, Eigen::Matrix3Xd positions, std::vector<HBondIDData> hbondAtomsMap, std::vector<AromaticIDData> aromaticAtomsMap);
    ~Context();

    // Copy constructor
    Context(const Context& other);
    // Copy assignment operator
    Context& operator=(const Context& other);

    void initializeContext();

    void computeContacts(size_t cellID);
    void computeBackboneAngles();
    void computeSidechainAngles(std::unordered_map<size_t, std::vector<TorsionIDData>> rotatingSCAtomsMap);
    void computeAromaticAngles();
    void computeHBondStates();
    void computeEnergy(Eigen::MatrixXd muPotential, std::vector<BBTorsionParamArray> bbTorsionParamVector, std::vector<SCTorsionParamArray> scTorsionParamVector, std::map<int, double> hbondPotential, HBondParamArray seq_based_hbondPotential,std::map<int, double> aromaticPotential);
    void removePosiotionsByAtomID(const std::vector<int>& atomIDs); // Not being used
    std::tuple<Eigen::Matrix3Xd, Eigen::Matrix3Xd> splitBackboneSidechainPositions(std::vector<size_t> backboneAtomIDs, std::vector<size_t> sidechainAtomIDs);
    void updateContext(std::vector<size_t> movedAtomIDs, std::vector<size_t> movedResidueIDs, Eigen::Matrix3Xd new_positions);
    Eigen::MatrixXd getContacts() const;
    Eigen::MatrixXd getBinaryContacts(double cutoff) const;
    std::array<Eigen::VectorXi, 4> getAngleIndices();
    std::vector<std::array<double, 4>> getSideChainTorsion();
    Eigen::Matrix3Xd getPositions() const;
    CellList getCellList() const;
    std::vector<size_t> setMovedAtomIDs();

    // Angles
    std::vector<double> getPhiAngles(std::vector<size_t> atomIDs);
    std::vector<double> getPsiAngles(std::vector<size_t> atomIDs);
    std::vector<double> getPlaneAngles(std::vector<size_t> atomIDs);
    std::vector<double> getBisectAngles(std::vector<size_t> atomIDs);

    // Energy
    Energy getEnergy() const;
    
private:
    Eigen::MatrixXd contacts; // Member variable to store contacts - this contact map is based on the cell list data structure (not all pairwise distances between atoms are computed)
    Eigen::Matrix3Xd positions; // Store the positions of atoms
    std::vector<size_t> movedAtomIDs; // Store the atoms that have moved
    Eigen::VectorXd backbone_bondlengths; // Store the bond lengths of the backbone atoms
    Eigen::VectorXd phi_angles; // Store the phi angles of the backbone atoms
    Eigen::VectorXd psi_angles; // Store the psi angles of the backbone atoms
    Eigen::VectorXd plane_angles; // Store the plane angles of the backbone atoms
    Eigen::VectorXd bisect_angles; // Store the bisect angles of the backbone atoms
    std::vector<double> aromatic_angles; // Store the aromatic angles
    Eigen::VectorXi phi_indices; // Store the indices of the phi angles
    Eigen::VectorXi psi_indices; // Store the indices of the psi angles
    Eigen::VectorXi plane_indices; // Store the indices of the plane angles
    Eigen::VectorXi bisect_indices; // Store the indices of the bisect angles
    std::vector<std::array<double, 4>> sidechain_torsions; // Store the sidechain torsions
    std::vector<std::vector<HBondState>> hbondStatesVector; // Store the hydrogen bond states
    Energy energy; // Store the energy values
    
    Topology& topology; // Reference to Topology object
    int numAtoms; // Member variable to store the number of atoms
    int numResidues; // Member variable to store the number of residues
    CellList cellList; // Member variable to store the cell list
    std::vector<BackBoneAtoms> backboneAtomsVector; // Store the backbone atomIDs
    std::vector<DihedralAtoms> sidechainAtomsVector; // Store the sidechain atomIDs
    std::vector<HBondIDData> hbondAtomsMap ; // Store the hydrogen bond atomIDs
    std::vector<AromaticIDData> aromaticAtomsMap; // Store the aromatic atomIDs

    void updatePositions(std::vector<size_t> movedAtomIDs, Eigen::Matrix3Xd new_positions);
    void updateContacts();
    void updateHBondStates(std::vector<size_t> movedResidueIDs);
    void updateAromaticAngles();
    void updateContacts2();
    void updateContacts3();
    void updateAngles(std::vector<size_t> movedResidueIDs);

    // Angle utils
    /*
    double computeDihedralAngle(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c, const Eigen::Vector3d& d);
    double computePlaneAngle(const Eigen::Vector3d& n1, const Eigen::Vector3d& ca1, const Eigen::Vector3d& o1,
                         const Eigen::Vector3d& n2, const Eigen::Vector3d& ca2, const Eigen::Vector3d& o2);
    double computeBisectAngle(const Eigen::Vector3d& n1, const Eigen::Vector3d& ca1, const Eigen::Vector3d& o1,
                         const Eigen::Vector3d& n2, const Eigen::Vector3d& ca2, const Eigen::Vector3d& o2);
    */
    Eigen::VectorXi convertToIndices(double cutoff, Eigen::VectorXd angles) const;
};

#endif // CONTEXT_H