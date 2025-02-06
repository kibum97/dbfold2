#include "system.h"
#include <iostream>
#include <string>
#include <algorithm>
#include "utils/hbond.h"

System::System(Topology& topology, MCPUForceField& forcefield)
    : topology(topology), forcefield(forcefield), muParamMatrix(0,0), smogTypeVector(0), numAtoms(0) {
    // Constructor implementation
    numAtoms = topology.getNumAtoms();
    numResidues = topology.getNumResidues();
    mupotential_parameters();
    torsion_parameters();
}

System::~System() {
    // Destructor implementation
}

// Topology Related Functions
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

std::unordered_map<size_t, std::vector<size_t>> System::getRotatingAtomsMap(Topology& topology) {
    // Get the map which stores value of atoms that move with respect to the reference atom
    std::unordered_map<size_t, std::vector<size_t>> rotatingAtomsMap;
    for (const auto& refAtom : topology.atoms) {
        if (refAtom.atomName == "N" || refAtom.atomName == "CA" || refAtom.atomName == "C") {
            size_t refResidueID = refAtom.resID;
            for (const auto& candidateAtom : topology.atoms) {
                if (candidateAtom.resID > refResidueID) {
                    rotatingAtomsMap[refAtom.atomID].push_back(candidateAtom.atomID);
                } else if(candidateAtom.resID == refResidueID) {
                    if (refAtom.atomName == "N") {
                        if (candidateAtom.atomName != "N" && candidateAtom.atomName != "H") {
                            rotatingAtomsMap[refAtom.atomID].push_back(candidateAtom.atomID);
                        }
                    } else if (refAtom.atomName == "CA") {
                        if (candidateAtom.atomName == "C" || candidateAtom.atomName == "O" || candidateAtom.atomName == "OXT") {
                            rotatingAtomsMap[refAtom.atomID].push_back(candidateAtom.atomID);
                        }
                    }                  
                }
            }
        }
    }
    return rotatingAtomsMap;
}

std::unordered_map<size_t, std::vector<TorsionIDData>> System::getSideChainRotatingAtomsMap(Topology& topology, std::unordered_map<std::string, std::vector<TorsionData>> torsion_map) {
    // Get the map which stores value of atoms that move with respect to the reference atom
    std::unordered_map<size_t, std::vector<TorsionIDData>> rotatingSCAtomsMap;
    for (const auto& refResidue: topology.residues) {
        std::string resName = refResidue.resName;
        size_t resID = refResidue.resID;
        TorsionIDData torsionidData;
        if (torsion_map.find(resName) == torsion_map.end()) {
            continue;
        }
        for (const auto& torsionData: torsion_map.at(resName)) {
            torsionidData.torsion_id = torsionData.torsion_id;
            if (torsionData.torsion_atoms.size() == 0) {
                continue;
            }
            size_t torsion_atom_index = 0;
            for (const auto& atomName: torsionData.torsion_atoms) {
                for (const auto& atom: refResidue.atoms) {
                    if (atom.atomName == atomName) {
                        torsionidData.torsion_atomIDs[torsion_atom_index] = atom.atomID;
                        torsion_atom_index++;
                    }
                }
            }
            for (const auto& atomName: torsionData.rotating_atoms) {
                for (const auto& atom: refResidue.atoms) {
                    if (atom.atomName == atomName) {
                        torsionidData.rotating_atomIDs.push_back(atom.atomID);
                    }
                }
            }
            rotatingSCAtomsMap[refResidue.resID].push_back(torsionidData);
        }
    }
    return rotatingSCAtomsMap;
}

std::vector<AromaticIDData> System::getAromaticAtomsMap(Topology& topology, std::unordered_map<std::string, AromaticData> aromatic_map) {
    // Get the map which stores value of atoms that move with respect to the reference atom
    std::vector<AromaticIDData> aromaticAtomsMap;
    size_t aromatid_id = 0;
    for (const auto& refResidue: topology.residues) {
        std::string resName = refResidue.resName;
        size_t resID = refResidue.resID;
        AromaticIDData aromaticIDData;
        if (aromatic_map.find(resName) == aromatic_map.end()) {
            continue;
        }
        size_t aromatic_atom_index = 0;
        for (const auto& atomName: aromatic_map.at(resName).aromatic_ring_atoms) {
            for (const auto& atom: refResidue.atoms) {
                if (atom.atomName == atomName) {
                    aromaticIDData.aromatic_atomIDs[aromatic_atom_index] = atom.atomID;
                    aromatic_atom_index++;
                    break;
                }
            }
        }
        aromaticIDData.aromatic_id = aromatid_id;
        aromaticAtomsMap.push_back(aromaticIDData);
        aromatid_id++;
    }
    return aromaticAtomsMap;
}

std::vector<HBondIDData> System::getHBondAtomsMap(Topology& topology, HBondData hbond_map) {
    // Get the map which stores value of atoms that move with respect to the reference atom
    std::vector<HBondIDData> hbondAtomsMap;
    size_t hbond_id = 0;
    for (size_t resID = 1; resID < topology.getNumResidues() - 1; ++resID) {
        HBondIDData hbondIDData;
        for (size_t donor_atom_index = 0; donor_atom_index < 7; ++donor_atom_index) {
            std::pair<int, std::string> atomInfo = hbond_map.donor_atoms[donor_atom_index];
            int offset = atomInfo.first;
            std::string atomName = atomInfo.second;
            Residue refResidue = topology.residues[resID + offset];
            std::cout << "Residue ID: " << resID + offset << " Atom Name: " << atomName << std::endl;
            for (const auto& atom: refResidue.atoms) {
                if (atom.atomName == atomName) {
                    hbondIDData.donor_atomIDs[donor_atom_index] = atom.atomID;
                    break;
                }
            }
        }
        for (size_t acceptor_atom_index = 0; acceptor_atom_index < 8; ++acceptor_atom_index) {
            std::pair<int, std::string> atomInfo = hbond_map.acceptor_atoms[acceptor_atom_index];
            int offset = atomInfo.first;
            std::string atomName = atomInfo.second;
            Residue refResidue = topology.residues[resID + offset];
            for (const auto& atom: refResidue.atoms) {
                if (atom.atomName == atomName) {
                    hbondIDData.acceptor_atomIDs[acceptor_atom_index] = atom.atomID;
                    break;
                }
            }
        }
        hbondIDData.hbond_id = hbond_id;
        hbondAtomsMap.push_back(hbondIDData);
        hbond_id++;
    }
    return hbondAtomsMap;
}

std::unordered_map<size_t, std::vector<RotamerData>> System::getRotamerMap(Topology& topology, std::string rotamer_data_file) {
    // Get the map which stores value of atoms that move with respect to the reference atom
    std::unordered_map<size_t, std::vector<RotamerData>> rotamerIDMap;
    std::unordered_map<std::string, std::vector<RotamerData>> rotamerMap = readRotamerData(rotamer_data_file);
    for (const auto& refResidue: topology.residues) {
        size_t resID = refResidue.resID;
        std::string resName = refResidue.resName;
        for (const auto& rotamerData: rotamerMap[resName]) {
            rotamerIDMap[resID].push_back(rotamerData);
        }
    }
    return rotamerIDMap;
}

// Contact Energy Related Functions
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
            //int index = computeHBondindex(0, 0, 0, 0, 0, 0, 0);
            int index = 0;
            // Compute hydrogen bond potential
            double value = hbondPotential[index];
        }
    }
}

// Backbone Torsion related functions
void System::torsion_parameters() {
    // Compute torsional parameters implementation
    BBTorsionArray bbTorsions = forcefield.getBackBoneTorsions();
    SCTorsionArray scTorsions = forcefield.getSideChainTorsions();
    std::vector<int> resTypeVector;
    for (int i = 0; i < numResidues; ++i) {
        int aa_index = lookupAminoAcidIndex(topology.residues[i].resName);
        resTypeVector.push_back(aa_index);
    }
    for (int i = 0; i < numResidues - 2; ++i) {
            BBTorsionParamArray bbTorsionParams = bbTorsions[resTypeVector[i]][resTypeVector[i+1]][resTypeVector[i+2]];
            SCTorsionParamArray scTorsionParams = scTorsions[resTypeVector[i]][resTypeVector[i+1]][resTypeVector[i+2]];
            bbTorsionParamVector.push_back(bbTorsionParams);
            scTorsionParamVector.push_back(scTorsionParams);
    }
}

int System::lookupAminoAcidIndex(const std::string& residueName) {
    // Define the map within the function
    static const std::unordered_map<std::string, int> aminoAcidIndex = {
        {"ALA", 0}, {"ARG", 1}, {"ASN", 2}, {"ASP", 3}, {"CYS", 4},
        {"GLN", 5}, {"GLU", 6}, {"GLY", 7}, {"HIS", 8}, {"ILE", 9},
        {"LEU", 10}, {"LYS", 11}, {"MET", 12}, {"PHE", 13}, {"PRO", 14},
        {"SER", 15}, {"THR", 16}, {"TRP", 17}, {"TYR", 18}, {"VAL", 19}
    };    
    // Convert input to uppercase
    std::string key = residueName;
    std::transform(key.begin(), key.end(), key.begin(), ::toupper);
    // Look up in the map
    auto it = aminoAcidIndex.find(key);
    if (it != aminoAcidIndex.end()) {
        return it->second; // Return the value if found
    } else {
        std::cerr << "Error: Invalid amino acid key '" << key << "'\n";
        return -1; // Sentinel value indicating failure
    }
}


// Accessors
Eigen::MatrixXd System::getMuParameters() const {
    return muParamMatrix;
}

std::vector <int> System::getSmogType() const {
    return smogTypeVector;
}

std::vector<BBTorsionParamArray> System::getBBTorsionParameters() const {
    return bbTorsionParamVector;
}

std::vector<SCTorsionParamArray> System::getSCTorsionParameters() const {
    return scTorsionParamVector;
}