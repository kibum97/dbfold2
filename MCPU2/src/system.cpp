#include "system.h"
#include <iostream>

System::System(Topology& topology, MCPUForceField& forcefield)
    : topology(topology), forcefield(forcefield), muParamMatrix(0,0), smogTypeVector(0), numAtoms(0) {
    // Constructor implementation
    numAtoms = topology.getNumAtoms();
    mupotential_parameters();
}

System::~System() {
    // Destructor implementation
}

std::tuple<std::vector<size_t>, std::vector<size_t>> System::splitBackboneSidechain(Topology& topology) {
    // Split the backbone and side chain atoms
    std::vector<size_t> backboneAtomIDs;
    std::vector<size_t> sidechainAtomIDs;
    for (const auto& atom : topology.atoms) {
        if (atom.atomName == "N" || atom.atomName == "CA" || atom.atomName == "C") {
            backboneAtomIDs.push_back(atom.atomID);
        } else {
            sidechainAtomIDs.push_back(atom.atomID);
        }
    }
    return std::make_tuple(backboneAtomIDs, sidechainAtomIDs);
}

void System::mupotential_parameters() {
    // Compute mu potential parameters implementation
    auto [smogTypeVector, remove_atom_ids] = forcefield.getSmogType(topology);
    Eigen::MatrixXd muPotential = forcefield.getMuPotential();
    //std::unordered_map <std::string, int> smogTypeMap = forcefield.getSmogTypeMap();
    //Remove atoms with no smog type
    if (remove_atom_ids.size() > 0) {
        topology.removeAtomsByID(remove_atom_ids);
        numAtoms = topology.getNumAtoms();
        std::cout << "Number of atoms after removing atoms with no smog type: " << numAtoms << std::endl;   
    }
    muParamMatrix = Eigen::MatrixXd::Zero(numAtoms, numAtoms);
    for (int i = 0; i < numAtoms; ++i) {
        for (int j = i + 1; j < numAtoms; ++j) {
            muParamMatrix(i, j) = muPotential(smogTypeVector[i], smogTypeVector[j]);
        }
    }
}

void System::hbond_parameters() {
    // Compute hydrogen bond parameters implementation
    std::map<int, double> hbondPotential = forcefield.getHBondPotential();
    HBondConfig hbondConfig = forcefield.getHBondConfig();

    for (int i = 1; i < numAtoms - 1; ++i) {
        for (int j = i + 1; j < numAtoms - 1; ++j) {
            // Compute residue index that contributes to hydrogen bond donor
            int index = forcefield.computeHBondindex(0, 0, 0, 0, 0, 0, 0);
            // Compute hydrogen bond potential
            double value = hbondPotential[index];
        }
    }
}



Eigen::MatrixXd System::getMuParameters() const {
    return muParamMatrix;
}

std::vector <int> System::getSmogType() const {
    return smogTypeVector;
}