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

void Topology::removeAtomsByID(const std::vector<int>& atomIDs) {
    // Remove from atoms vector
    atoms.erase(std::remove_if(atoms.begin(), atoms.end(),
        [&atomIDs](const Atom& atom) { return std::find(atomIDs.begin(), atomIDs.end(), atom.atomID) != atomIDs.end(); }), atoms.end());

    // Remove from residues vector
    for (auto& residue : residues) {
        residue.atoms.erase(std::remove_if(residue.atoms.begin(), residue.atoms.end(),
            [&atomIDs](const Atom& atom) { return std::find(atomIDs.begin(), atomIDs.end(), atom.atomID) != atomIDs.end(); }), residue.atoms.end());
    }

    // Remove empty residues and chains
    residues.erase(std::remove_if(residues.begin(), residues.end(),
        [](const Residue& residue) { return residue.atoms.empty(); }), residues.end());
    for (auto& chain : chains) {
        chain.residues.erase(std::remove_if(chain.residues.begin(), chain.residues.end(),
            [](const Residue& residue) { return residue.atoms.empty(); }), chain.residues.end());
    }

    // Update topology counts
    updateTopology();
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