#include "topology.h"
#include <iostream>

int Topology::getNumAtoms() const {
    return numAtoms;
}

int Topology::getNumResidues() const {
    return numResidues;
}

int Topology::getNumChains() const {
    return numChains;
}

void Topology::setPeriodicBox(const Eigen::Matrix3d& periodicBox) {
    this->periodicBox = periodicBox;
}

Eigen::Matrix3d Topology::getPeriodicBox() const {
    return this->periodicBox;
}

void Topology::updateTopology() {
    numAtoms = atoms.size();
    numResidues = residues.size();
    numChains = chains.size();
}

void Topology::printTopology() {
    std::cout << "Number of atoms: " << numAtoms << std::endl;
    std::cout << "Number of residues: " << numResidues << std::endl;
    std::cout << "Number of chains: " << numChains << std::endl;
    std::cout << "Periodic box size: " << periodicBox << std::endl;
}