#ifndef CONTEXT_H
#define CONTEXT_H

#include <vector>
#include <Eigen/Dense>
#include <unordered_set>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "topology.h"
#include "cell_list.h"

class Context {
public:
    Context(Topology& topology, Eigen::Matrix3Xd positions);
    ~Context();

    void computeContacts(size_t cellID);
    void removePosiotionsByAtomID(const std::vector<int>& atomIDs); // Not being used
    std::tuple<Eigen::Matrix3Xd, Eigen::Matrix3Xd> splitBackboneSidechainPositions(std::vector<size_t> backboneAtomIDs, std::vector<size_t> sidechainAtomIDs);
    void updateContext(std::vector<size_t> movedAtomIDs, Eigen::Matrix3Xd new_positions);
    Eigen::MatrixXd getContacts() const;
    Eigen::MatrixXd getBinaryContacts(double cutoff) const;
    Eigen::Matrix3Xd getPositions() const;
    CellList getCellList() const;
    std::vector<size_t> setMovedAtomIDs();
    
private:
    Eigen::MatrixXd contacts; // Member variable to store contacts - this contact map is based on the cell list data structure (not all pairwise distances between atoms are computed)
    Eigen::Matrix3Xd positions; // Store the positions of atoms
    std::vector<size_t> movedAtomIDs; // Store the atoms that have moved
    Eigen::VectorXd backbone_bondlengths; // Store the bond lengths of the backbone atoms
    Eigen::VectorXd backbone_diheral_angles; // Store the dihedral angles of the backbone atoms
    Topology& topology; // Reference to Topology object
    int numAtoms; // Member variable to store the number of atoms
    CellList cellList; // Member variable to store the cell list

    void updatePositions(std::vector<size_t> movedAtomIDs, Eigen::Matrix3Xd new_positions);
    void updateContacts();
};

#endif // CONTEXT_H