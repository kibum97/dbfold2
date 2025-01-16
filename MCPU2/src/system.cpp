#include "system.h"
#include <iostream>

System::System(Topology& topology, MCPUForceField& forcefield)
    : topology(topology), forcefield(forcefield), muParamMatrix(0,0), smogTypeVector(0), numAtoms(0) {
    // Constructor implementation
}

System::~System() {
    // Destructor implementation
}

void System::initialize() {
    // Initialize context implementation
    numAtoms = topology.getNumAtoms();
    muParamMatrix = Eigen::MatrixXd::Zero(numAtoms, numAtoms);
    smogTypeVector.resize(numAtoms);
    mupotential_parameters();
}

void System::mupotential_parameters() {
    // Compute mu potential parameters implementation
    std::string sidechain_key, backbone_key;
    Eigen::MatrixXd muPotential = forcefield.getMuPotential();
    std::unordered_map <std::string, int> smogTypeMap = forcefield.getSmogTypeMap();
    // Assign Smog Type
    for (int i = 0; i < numAtoms; ++i) {
        sidechain_key = topology.getAtom(i).resName + std::string("_") + topology.getAtom(i).atomName;
        backbone_key = std::string("XXX_") + topology.getAtom(i).atomName;
        if (smogTypeMap.find(sidechain_key) != smogTypeMap.end()) {
            int smogType = smogTypeMap[sidechain_key];
            smogTypeVector[i] = smogType;
        } else if (smogTypeMap.find(backbone_key) != smogTypeMap.end()) {
            int smogType = smogTypeMap[backbone_key];
            smogTypeVector[i] = smogType;
        } else {
            std::cerr << "Error: Smog Type not found for atom " << topology.getAtom(i).atomID << std::endl;
        }
    }

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