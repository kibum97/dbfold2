#include "context.h"
#include "utils/geometry.h"
#include "utils/hbond.h"

//TODO: Need to specify way to define cell size
//TODO: vector of angles - need to think about the length of each vector
Context::Context(Topology& topology, Eigen::Matrix3Xd positions, std::vector<HBondIDData> hbondAtomsMap, std::vector<AromaticIDData> aromaticAtomsMap)
    : topology(topology), numAtoms(topology.getNumAtoms()), numResidues(topology.getNumResidues()), positions(positions), cellList(positions, 5.0, topology.getPeriodicBox()), contacts(0, 0) {
    // Store value maps to local variables
    this->hbondAtomsMap = hbondAtomsMap;
    this->aromaticAtomsMap = aromaticAtomsMap;    
    // Constructor implementation
    contacts = Eigen::MatrixXd::Constant(numAtoms, numAtoms, 10.0);
    phi_angles = Eigen::VectorXd::Zero(numResidues-2);
    psi_angles = Eigen::VectorXd::Zero(numResidues-2);
    plane_angles = Eigen::VectorXd::Zero(numResidues-2);
    bisect_angles = Eigen::VectorXd::Zero(numResidues-2);
    phi_indices = Eigen::VectorXi::Zero(numResidues-2);
    psi_indices = Eigen::VectorXi::Zero(numResidues-2);
    plane_indices = Eigen::VectorXi::Zero(numResidues-2);
    bisect_indices = Eigen::VectorXi::Zero(numResidues-2);
}

Context::~Context() {
    // Destructor implementation
}

// Implement the copy constructor
Context::Context(const Context& other)
    : topology(other.topology),
      positions(other.positions),
      hbondAtomsMap(other.hbondAtomsMap),
      aromaticAtomsMap(other.aromaticAtomsMap),
      cellList(other.cellList) {
    // Perform deep copy of other resources if necessary
    contacts = other.contacts;
    phi_angles = other.phi_angles;
    psi_angles = other.psi_angles;
    plane_angles = other.plane_angles;
    bisect_angles = other.bisect_angles;
    phi_indices = other.phi_indices;
    psi_indices = other.psi_indices;
    plane_indices = other.plane_indices;
    bisect_indices = other.bisect_indices;
    sidechain_torsions = other.sidechain_torsions;
    hbondStatesVector = other.hbondStatesVector;
    aromatic_angles = other.aromatic_angles;
    energy = other.energy;
    backboneAtomsVector = other.backboneAtomsVector;
    sidechainAtomsVector = other.sidechainAtomsVector;
    numAtoms = other.numAtoms;
    numResidues = other.numResidues;
}
// Implement the copy assignment operator
Context& Context::operator=(const Context& other) {
    if (this == &other) {
        return *this; // Handle self-assignment
    }
    topology = other.topology;
    positions = other.positions;
    hbondAtomsMap = other.hbondAtomsMap;
    aromaticAtomsMap = other.aromaticAtomsMap;
    cellList = other.cellList;
    // Perform deep copy of other resources if necessary
    contacts = other.contacts;
    phi_angles = other.phi_angles;
    psi_angles = other.psi_angles;
    plane_angles = other.plane_angles;
    bisect_angles = other.bisect_angles;
    phi_indices = other.phi_indices;
    psi_indices = other.psi_indices;
    plane_indices = other.plane_indices;
    bisect_indices = other.bisect_indices;
    sidechain_torsions = other.sidechain_torsions;
    hbondStatesVector = other.hbondStatesVector;
    aromatic_angles = other.aromatic_angles;
    energy = other.energy;
    backboneAtomsVector = other.backboneAtomsVector;
    sidechainAtomsVector = other.sidechainAtomsVector;
    numAtoms = other.numAtoms;
    numResidues = other.numResidues;
    return *this;
}

void Context::initializeContext() {
    // Initialize context implementation
    // Compute contacts
    for (int i = 0; i < cellList.cells.size(); ++i) {
        computeContacts(i);
    }
    // Compute dihedral angles
    computeBackboneAngles();
    // Compute hydrogen bond states

}

void Context::computeContacts(size_t cellID) {
    // Compute contacts using cell list data type
    Cell cell = cellList.cells[cellID];
    for (auto atomID1 : cell.atomIDs) {
        for (int neighborID : cell.neighborIds) {
            Cell neighbor = cellList.cells[neighborID];
            for (auto atomID2 : neighbor.atomIDs) {
                double distance = (positions.col(atomID1) - positions.col(atomID2)).norm();
                contacts(atomID1, atomID2) = distance;
                contacts(atomID2, atomID1) = distance;
            }
        }
    }
}

void Context::computeEnergy(Eigen::MatrixXd muPotential, std::vector<BBTorsionParamArray> bbTorsionParamVector, std::vector<SCTorsionParamArray> scTorsionParamVector, std::map<int, double> hbondPotential, HBondParamArray seq_based_hbondPotential,std::map<int, double> aromaticPotential) {
    // Compute Mu potential energy
    auto mu_start = std::chrono::high_resolution_clock::now();
    Eigen::MatrixXd binaryContacts = getBinaryContacts(5.0);
    Eigen::MatrixXd result = muPotential.array() * binaryContacts.array();
    energy.mu_energy = result.sum();
    std::cout << "Mu potential energy: " << energy.mu_energy << std::endl;
    auto mu_end = std::chrono::high_resolution_clock::now();

    // Compute Backbone Torsion energy
    double torsion_energy = 0.0;
    for (int resid = 0; resid < topology.getNumResidues()-2; ++resid) {
        int phi_index = phi_indices[resid];
        int psi_index = psi_indices[resid];
        int plane_index = plane_indices[resid];
        int bisect_index = bisect_indices[resid];
        torsion_energy += bbTorsionParamVector[resid][phi_index][psi_index][plane_index][bisect_index];
    }
    torsion_energy = torsion_energy/1000.; 
    energy.bb_torsion_energy = torsion_energy;
    std::cout << "Backbone torsion energy: " << torsion_energy << std::endl;
    auto bb_torsion_end = std::chrono::high_resolution_clock::now();

    // Compute Sidechain torsion energy
    double sc_torsion_energy = 0.0;
    for (int i = 0; i < sidechain_torsions.size(); ++i) {
        std::array<int, 4> chi_indices = {0, 0, 0, 0};
        for (int j = 0; j < sidechain_torsions[i].size(); ++j) {
            double angle = sidechain_torsions[i][j];
            chi_indices[j] = static_cast<int>((normalizeAngle(angle) + 15.)/30.);
        }
        sc_torsion_energy += scTorsionParamVector[i][chi_indices[0]][chi_indices[1]][chi_indices[2]][chi_indices[3]];
    }
    sc_torsion_energy = sc_torsion_energy/1000.;
    energy.sc_torsion_energy = sc_torsion_energy;
    std::cout << "Sidechain torsion energy: " << sc_torsion_energy << std::endl;
    auto sc_torsion_end = std::chrono::high_resolution_clock::now();

    // Compute Hydrogen bond energy
    double hbond_energy = 0.0;
    // HBond - sequence specific
    for (auto hbond : hbondStatesVector) {
        for (auto state : hbond) {
            if (state.hbond_state) {
                // Sequence based
                hbond_energy += seq_based_hbondPotential[state.secondary_structure_label][state.donor_resType][state.acceptor_resType];
                // Geometry based
                int index = computeHBondindexfromAngles(state.secondary_structure_label, state.bisect_angles, state.plane_angles);
                hbond_energy += hbondPotential[index];
            }
        }
    }
    hbond_energy = hbond_energy/1000. * 2.;
    std::cout << "HBond energy: " << hbond_energy << std::endl;
    energy.hbond_energy = hbond_energy;
    auto hbond_end = std::chrono::high_resolution_clock::now();

    // Compute Aromatic energy
    double aromatic_energy = 0.0;
    for (auto angle : aromatic_angles) {
        if (angle > 89.9) {
            angle = 89.9;
        }
        int aromatic_index = static_cast<int>(angle/10.);
        aromatic_energy += aromaticPotential[aromatic_index];
    }
    aromatic_energy = aromatic_energy/1000.;
    energy.aromatic_energy = aromatic_energy;
    std::cout << "Aromatic energy: " << aromatic_energy << std::endl;
    auto aromatic_end = std::chrono::high_resolution_clock::now();

    // Compute total energy
    energy.total_energy = energy.mu_energy + energy.bb_torsion_energy + energy.sc_torsion_energy + hbond_energy + aromatic_energy;

    // Calculate elapsed time for each part
    std::chrono::duration<double> mu_elapsed = mu_end - mu_start;
    std::chrono::duration<double> bb_torsion_elapsed = bb_torsion_end - mu_end;
    std::chrono::duration<double> sc_torsion_elapsed = sc_torsion_end - bb_torsion_end;
    std::chrono::duration<double> hbond_elapsed = hbond_end - sc_torsion_end;
    std::chrono::duration<double> aromatic_elapsed = aromatic_end - hbond_end;

    // Print elapsed time for each part
    std::cout << "Elapsed time for Mu potential: " << mu_elapsed.count() << " seconds" << std::endl;
    std::cout << "Elapsed time for Backbone Torsion: " << bb_torsion_elapsed.count() << " seconds" << std::endl;
    std::cout << "Elapsed time for Sidechain Torsion: " << sc_torsion_elapsed.count() << " seconds" << std::endl;
    std::cout << "Elapsed time for HBond: " << hbond_elapsed.count() << " seconds" << std::endl;
    std::cout << "Elapsed time for Aromatic: " << aromatic_elapsed.count() << " seconds" << std::endl;

}

void Context::computeBackboneAngles() {
    // Lookup backbone atoms
    for (auto residue: topology.residues) {
        BackBoneAtoms backboneAtoms = {};
        for (const auto& atom : residue.atoms) {
            if (atom.atomName == "N") {
                backboneAtoms[0] = atom.atomID;
            } else if (atom.atomName == "CA") {
                backboneAtoms[1] = atom.atomID;
            } else if (atom.atomName == "C") {
                backboneAtoms[2] = atom.atomID;
            } else if (atom.atomName == "O") {
                backboneAtoms[3] = atom.atomID;
            }
        }
        backboneAtomsVector.push_back(backboneAtoms);
    }
    std::cout << "Backbone atoms saved" << std::endl;
    // Compute phi angles
    std::cout << "Computing phi angles" << std::endl;
    for (int i = 1; i < numResidues - 1; ++i) {
        double phi = computeDihedralAngle(positions.col(backboneAtomsVector[i-1][2]), positions.col(backboneAtomsVector[i][0]), positions.col(backboneAtomsVector[i][1]), positions.col(backboneAtomsVector[i][2]));
        phi_angles[i-1] = phi * (180.0 / M_PI);
        phi_indices[i-1] = phi/60.0;
    }
    std::cout << "Phi angles: " << phi_angles << std::endl;
    // Compute psi angles
    std::cout << "Computing psi angles" << std::endl;
    for (int i = 1; i < numResidues - 1; ++i) {
        double psi = computeDihedralAngle(positions.col(backboneAtomsVector[i][0]), positions.col(backboneAtomsVector[i][1]), positions.col(backboneAtomsVector[i][2]), positions.col(backboneAtomsVector[i+1][0]));
        psi_angles[i-1] = psi * (180.0 / M_PI);
        psi_indices[i-1] = psi/60.0;
    }
    std::cout << "Psi angles: " << psi_angles << std::endl;
    // Compute plane angles
    std::cout << "Computing plane angles" << std::endl;
    for (int i = 0; i < numResidues - 2; ++i) {
        double plane = computePlaneAngle(positions.col(backboneAtomsVector[i][0]), positions.col(backboneAtomsVector[i][1]), positions.col(backboneAtomsVector[i][3]),
                                         positions.col(backboneAtomsVector[i+2][0]), positions.col(backboneAtomsVector[i+2][1]), positions.col(backboneAtomsVector[i+2][3]));
        plane_angles[i] = plane * (180.0 / M_PI);
        plane_indices[i] = plane/30.0;
    }
    std::cout << "Plane angles: " << plane_angles << std::endl;
    // Compute bisect angles
    std::cout << "Computing bisect angles" << std::endl;
    for (int i = 0; i < numResidues - 2; ++i) {
        double bisect = computeBisectAngle(positions.col(backboneAtomsVector[i][0]), positions.col(backboneAtomsVector[i][1]), positions.col(backboneAtomsVector[i][3]),
                                           positions.col(backboneAtomsVector[i+2][0]), positions.col(backboneAtomsVector[i+2][1]), positions.col(backboneAtomsVector[i+2][3]));
        bisect_angles[i] = bisect * (180.0 / M_PI);
        bisect_indices[i] = bisect/30.0;
    }
    std::cout << "Bisect angles: " << bisect_angles << std::endl;
}

void Context::computeSidechainAngles(std::unordered_map<size_t, std::vector<TorsionIDData>> rotatingSCAtomsMap) {
    // Compute side chain angles
    for (auto scAtomsMap : rotatingSCAtomsMap) {
        size_t resID = scAtomsMap.first;
        std::vector<TorsionIDData> torsionidData = scAtomsMap.second;
        std::array<double, 4> chi_angles = {0., 0., 0., 0.};
        for (auto data: torsionidData) {
            double chi_angle;
            std::array<size_t, 4> atomIDs = data.torsion_atomIDs;
            chi_angles[data.torsion_id] = computeDihedralAngle(positions.col(atomIDs[0]), positions.col(atomIDs[1]), positions.col(atomIDs[2]), positions.col(atomIDs[3]));
        }
        sidechain_torsions.push_back(chi_angles);
    }
}

void Context::computeAromaticAngles() {
    aromatic_angles.clear();
    // Compute aromatic angles
    for (auto i = 0; i < aromaticAtomsMap.size(); ++i) {
        AromaticIDData aromaticIDData1 = aromaticAtomsMap[i];
        std::array<size_t, 3> atomIDs1 = aromaticIDData1.aromatic_atomIDs;
        Eigen::Vector3d center1 = computeCenter({positions.col(atomIDs1[0]), positions.col(atomIDs1[1]), positions.col(atomIDs1[2])});
        for (auto j = i + 1; j < aromaticAtomsMap.size(); ++j) {
            AromaticIDData aromaticIDData2 = aromaticAtomsMap[j];
            std::array<size_t, 3> atomIDs2 = aromaticIDData2.aromatic_atomIDs;
            Eigen::Vector3d center2 = computeCenter({positions.col(atomIDs2[0]), positions.col(atomIDs2[1]), positions.col(atomIDs2[2])});
            if ((center1 - center2).norm() < 7.0) {
                double angle = computePlaneAngle(positions.col(atomIDs1[0]), positions.col(atomIDs1[1]), positions.col(atomIDs1[2]),
                                             positions.col(atomIDs2[0]), positions.col(atomIDs2[1]), positions.col(atomIDs2[2]));
                aromatic_angles.push_back(angle);
            }
        }
    }
}

void Context::computeHBondStates() {
    // Compute hydrogen bond states
    hbondStatesVector.resize(hbondAtomsMap.size());
    for (auto i = 0; i < hbondAtomsMap.size(); ++i) {
        HBondIDData hbondIDData_donor = hbondAtomsMap[i];
        std::array<size_t, 7> donor_atomIDs = hbondIDData_donor.donor_atomIDs;
        size_t donorResID = topology.atoms[donor_atomIDs[0]].resID;
        std::cout << "Donor ResID: " << donorResID << " index: " << i << std::endl;
        for (auto j = 0; j < hbondAtomsMap.size(); ++j) {
            HBondState hbondState;
            HBondIDData hbondIDData_acceptor = hbondAtomsMap[j];
            std::array<size_t, 8> acceptor_atomIDs = hbondIDData_acceptor.acceptor_atomIDs;
            size_t acceptorResID = topology.atoms[acceptor_atomIDs[0]].resID;
            double bisect_angle1 = computeBisectAngle(positions.col(donor_atomIDs[0]), positions.col(donor_atomIDs[1]), positions.col(donor_atomIDs[3]),
                                                      positions.col(acceptor_atomIDs[4]), positions.col(acceptor_atomIDs[2]), positions.col(acceptor_atomIDs[1]));
            double plane_angle1 = computePlaneAngle(positions.col(donor_atomIDs[0]), positions.col(donor_atomIDs[1]), positions.col(donor_atomIDs[3]),
                                                    positions.col(acceptor_atomIDs[4]), positions.col(acceptor_atomIDs[2]), positions.col(acceptor_atomIDs[1]));
            double bisect_angle2 = computeBisectAngle(positions.col(donor_atomIDs[5]), positions.col(donor_atomIDs[4]), positions.col(donor_atomIDs[2]),
                                                      positions.col(acceptor_atomIDs[3]), positions.col(acceptor_atomIDs[5]), positions.col(acceptor_atomIDs[6]));
            double plane_angle2 = computePlaneAngle(positions.col(donor_atomIDs[5]), positions.col(donor_atomIDs[4]), positions.col(donor_atomIDs[2]),
                                                    positions.col(acceptor_atomIDs[3]), positions.col(acceptor_atomIDs[5]), positions.col(acceptor_atomIDs[6]));
            double bisect_angle3 = computeBisectAngle(positions.col(donor_atomIDs[2]), positions.col(donor_atomIDs[0]), positions.col(donor_atomIDs[1]),
                                                      positions.col(acceptor_atomIDs[2]), positions.col(acceptor_atomIDs[1]), positions.col(acceptor_atomIDs[3]));
            double plane_angle3 = computePlaneAngle(positions.col(donor_atomIDs[2]), positions.col(donor_atomIDs[0]), positions.col(donor_atomIDs[1]),
                                                    positions.col(acceptor_atomIDs[2]), positions.col(acceptor_atomIDs[1]), positions.col(acceptor_atomIDs[3]));
            hbondState.donor_resType = topology.residues[donorResID].resType;
            hbondState.acceptor_resType = topology.residues[acceptorResID].resType;
            hbondState.bisect_angles = {bisect_angle1, bisect_angle2, bisect_angle3};
            hbondState.plane_angles = {plane_angle1, plane_angle2, plane_angle3};
            hbondState.hbond_state = checkHydrogenBond(donorResID, donor_atomIDs, acceptorResID, acceptor_atomIDs, positions);
            hbondState.secondary_structure_label = computeSecondaryStructure(donorResID, donor_atomIDs, acceptorResID, acceptor_atomIDs, positions);
            hbondStatesVector[i].push_back(hbondState);
        }
    }
}

void Context::removePosiotionsByAtomID(const std::vector<int>& atomIDs) {
    // Create a mask to keep track of columns to keep
    std::vector<bool> keep(positions.cols(), true);
    for (int index : atomIDs) {
        if (index >= 0 && index < positions.cols()) {
            keep[index] = false;
        }
    }
    // Count the number of columns to keep
    int newCols = std::count(keep.begin(), keep.end(), true);
    // Create a new matrix with the appropriate number of columns
    Eigen::Matrix3Xd newPositions(3, newCols);
    // Copy the columns that are to be kept
    for (int i = 0, j = 0; i < positions.cols(); ++i) {
        if (keep[i]) {
            newPositions.col(j++) = positions.col(i);
        }
    }
    // Replace the original matrix with the new matrix
    positions = newPositions;
}

std::tuple<Eigen::Matrix3Xd, Eigen::Matrix3Xd> Context::splitBackboneSidechainPositions(std::vector<size_t> backboneAtomIDs, std::vector<size_t> sidechainAtomIDs) {
    // Split the positions of backbone and side chain atoms
    Eigen::Matrix3Xd backbonePositions(3, backboneAtomIDs.size());
    Eigen::Matrix3Xd sidechainPositions(3, sidechainAtomIDs.size());
    for (int i = 0; i < backboneAtomIDs.size(); ++i) {
        backbonePositions.col(i) = positions.col(backboneAtomIDs[i]);
    }
    for (int i = 0; i < sidechainAtomIDs.size(); ++i) {
        sidechainPositions.col(i) = positions.col(sidechainAtomIDs[i]);
    }
    return std::make_tuple(backbonePositions, sidechainPositions);
}

void Context::updateContext(std::vector<size_t> movedAtomIDs, std::vector<size_t> movedResidueIDs, Eigen::Matrix3Xd new_positions) {
    auto context_start = std::chrono::high_resolution_clock::now();
    updatePositions(movedAtomIDs, new_positions);
    auto context_position = std::chrono::high_resolution_clock::now();
    updateAngles(movedResidueIDs);
    auto context_angle = std::chrono::high_resolution_clock::now();
    updateContacts2(); //TODO: Current bottleneck
    auto context_contact = std::chrono::high_resolution_clock::now();
    std::cout << "updating HBond states" << std::endl;
    updateHBondStates(movedResidueIDs);
    auto context_hbond = std::chrono::high_resolution_clock::now();
    std::cout << "updating Aromatic angles" << std::endl;
    updateAromaticAngles();
    auto context_aromatic = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> context_position_elapsed = context_position - context_start;
    std::chrono::duration<double> context_angle_elapsed = context_angle - context_position;
    std::chrono::duration<double> context_contact_elapsed = context_contact - context_angle;
    std::chrono::duration<double> context_hbond_elapsed = context_hbond - context_contact;
    std::chrono::duration<double> context_aromatic_elapsed = context_aromatic - context_hbond;

    // Print elapsed time for each part
    std::cout << "Elapsed time for updating positions: " << context_position_elapsed.count() << " seconds" << std::endl;
    std::cout << "Elapsed time for updating angles: " << context_angle_elapsed.count() << " seconds" << std::endl;
    std::cout << "Elapsed time for updating contacts: " << context_contact_elapsed.count() << " seconds" << std::endl;
    std::cout << "Elapsed time for updating HBond states: " << context_hbond_elapsed.count() << " seconds" << std::endl;
    std::cout << "Elapsed time for updating Aromatic angles: " << context_aromatic_elapsed.count() << " seconds" << std::endl;
}

void Context::updatePositions(std::vector<size_t> movedAtomIDs, Eigen::Matrix3Xd new_positions) {
    this->positions = new_positions;
    this->movedAtomIDs = movedAtomIDs;
}

void Context::updateContacts() {
    auto start_bottleneck = std::chrono::high_resolution_clock::now();
    std::unordered_set<size_t> movedCellIDs = cellList.updateAtoms(movedAtomIDs, positions); // TODO: Current bottleneck
    auto bn1_end = std::chrono::high_resolution_clock::now();
    for (auto cellID : movedCellIDs) {
        //if (cellID >= cellList.cells.size()) exit(1);
        computeContacts(cellID);
    }
    auto bn2_end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> bn1_elapsed = bn1_end - start_bottleneck;
    std::chrono::duration<double> bn2_elapsed = bn2_end - bn1_end;

    // Print elapsed time for each part
    std::cout << "Elapsed time for Bottleneck 1: " << bn1_elapsed.count() << " seconds" << std::endl;
    std::cout << "Elapsed time for Bottleneck 2: " << bn2_elapsed.count() << " seconds" << std::endl;
}


void Context::updateHBondStates(std::vector<size_t> movedResidueIDs) {
    for (auto donorResID : movedResidueIDs) {
        if (donorResID < 2 || donorResID > numResidues - 1) {
            std::cout << "Residue ID out of range" << std::endl;
            continue;
        }
        std::cout << "Size of hbondAtomsMap: " << hbondAtomsMap.size() << std::endl;
        HBondIDData hbondIDData_donor = hbondAtomsMap[donorResID-2];
        std::cout << "Donor ResID: " << donorResID << std::endl;
        std::array<size_t, 7> donor_atomIDs = hbondIDData_donor.donor_atomIDs;
        for (int i = 0; i < hbondAtomsMap.size(); ++i) {
            size_t acceptorResID = i + 2;
            HBondState new_hbondState;
            HBondIDData hbondIDData_acceptor = hbondAtomsMap[i];
            std::array<size_t, 8> acceptor_atomIDs = hbondIDData_acceptor.acceptor_atomIDs;
            double bisect_angle1 = computeBisectAngle(positions.col(donor_atomIDs[0]), positions.col(donor_atomIDs[1]), positions.col(donor_atomIDs[3]),
                                                      positions.col(acceptor_atomIDs[4]), positions.col(acceptor_atomIDs[2]), positions.col(acceptor_atomIDs[1]));
            double plane_angle1 = computePlaneAngle(positions.col(donor_atomIDs[0]), positions.col(donor_atomIDs[1]), positions.col(donor_atomIDs[3]),
                                                    positions.col(acceptor_atomIDs[4]), positions.col(acceptor_atomIDs[2]), positions.col(acceptor_atomIDs[1]));
            double bisect_angle2 = computeBisectAngle(positions.col(donor_atomIDs[5]), positions.col(donor_atomIDs[4]), positions.col(donor_atomIDs[2]),
                                                      positions.col(acceptor_atomIDs[3]), positions.col(acceptor_atomIDs[5]), positions.col(acceptor_atomIDs[6]));
            double plane_angle2 = computePlaneAngle(positions.col(donor_atomIDs[5]), positions.col(donor_atomIDs[4]), positions.col(donor_atomIDs[2]),
                                                    positions.col(acceptor_atomIDs[3]), positions.col(acceptor_atomIDs[5]), positions.col(acceptor_atomIDs[6]));
            double bisect_angle3 = computeBisectAngle(positions.col(donor_atomIDs[2]), positions.col(donor_atomIDs[0]), positions.col(donor_atomIDs[1]),
                                                      positions.col(acceptor_atomIDs[2]), positions.col(acceptor_atomIDs[1]), positions.col(acceptor_atomIDs[3]));
            double plane_angle3 = computePlaneAngle(positions.col(donor_atomIDs[2]), positions.col(donor_atomIDs[0]), positions.col(donor_atomIDs[1]),
                                                    positions.col(acceptor_atomIDs[2]), positions.col(acceptor_atomIDs[1]), positions.col(acceptor_atomIDs[3]));
            new_hbondState.bisect_angles = {bisect_angle1, bisect_angle2, bisect_angle3};
            new_hbondState.plane_angles = {plane_angle1, plane_angle2, plane_angle3};
            new_hbondState.hbond_state = checkHydrogenBond(donorResID, donor_atomIDs, acceptorResID, acceptor_atomIDs, positions);
            new_hbondState.secondary_structure_label = computeSecondaryStructure(donorResID, donor_atomIDs, acceptorResID, acceptor_atomIDs, positions);
            hbondStatesVector[donorResID - 2][i] = new_hbondState;
        }
    }
}

void Context::updateAromaticAngles() {
    this->computeAromaticAngles();
}

void Context::updateContacts2() {
    auto start_bottleneck = std::chrono::high_resolution_clock::now();
    auto bn1_end = std::chrono::high_resolution_clock::now();
    for (auto atomID1 : movedAtomIDs) {
        std::array<int, 3> cellIndex = cellList.getCellIndex(positions.col(atomID1));
        size_t cellID = cellList.flattenIndex(cellIndex);
        Cell cell = cellList.cells[cellID];
        for (auto neighborID : cell.neighborIds) {
            Cell neighbor = cellList.cells[neighborID];
            for (auto atomID2 : neighbor.atomIDs) {
                /*
                bool atom2_moved = std::find(movedAtomIDs.begin(), movedAtomIDs.end(), atomID2) != movedAtomIDs.end();
                if (!atom2_moved) {
                    double distance = (positions.col(atomID1) - positions.col(atomID2)).norm();
                    contacts(atomID1, atomID2) = distance;
                    contacts(atomID2, atomID1) = distance;
                }
                */
                double distance = (positions.col(atomID1) - positions.col(atomID2)).norm();
                contacts(atomID1, atomID2) = distance;
                contacts(atomID2, atomID1) = distance;
            }
        }
    }
    auto bn2_end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> bn1_elapsed = bn1_end - start_bottleneck;
    std::chrono::duration<double> bn2_elapsed = bn2_end - bn1_end;

    // Print elapsed time for each part
    std::cout << "Elapsed time for Bottleneck 1: " << bn1_elapsed.count() << " seconds" << std::endl;
    std::cout << "Elapsed time for Bottleneck 2: " << bn2_elapsed.count() << " seconds" << std::endl;
}

void Context::updateContacts3() {
    auto start_bottleneck = std::chrono::high_resolution_clock::now();
    auto bn1_end = std::chrono::high_resolution_clock::now();
    for (auto atomID1 : movedAtomIDs) {
        std::array<int, 3> cellIndex = cellList.getCellIndex(positions.col(atomID1));
        size_t cellID = cellList.flattenIndex(cellIndex);
        Cell cell = cellList.cells[cellID];
        for (auto neighborID : cell.neighborIds) {
            Cell neighbor = cellList.cells[neighborID];
            for (auto atomID2 : neighbor.atomIDs) {
                bool atom2_moved = std::find(movedAtomIDs.begin(), movedAtomIDs.end(), atomID2) != movedAtomIDs.end();
                if (!atom2_moved) {
                    double distance = (positions.col(atomID1) - positions.col(atomID2)).norm();
                    contacts(atomID1, atomID2) = distance;
                    contacts(atomID2, atomID1) = distance;
                }
                double distance = (positions.col(atomID1) - positions.col(atomID2)).norm();
                contacts(atomID1, atomID2) = distance;
                contacts(atomID2, atomID1) = distance;
            }
        }
    }
    auto bn2_end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> bn1_elapsed = bn1_end - start_bottleneck;
    std::chrono::duration<double> bn2_elapsed = bn2_end - bn1_end;

    // Print elapsed time for each part
    std::cout << "Elapsed time for Bottleneck 1: " << bn1_elapsed.count() << " seconds" << std::endl;
    std::cout << "Elapsed time for Bottleneck 2: " << bn2_elapsed.count() << " seconds" << std::endl;
}

void Context::updateAngles(std::vector<size_t> movedResidueIDs) {
    for (auto resID : movedResidueIDs) {
        BackBoneAtoms backboneAtoms = backboneAtomsVector[resID];
        // Compute phi angles
        if (resID > 0 && resID < numResidues - 1) {
            double phi = computeDihedralAngle(positions.col(backboneAtomsVector[resID-1][2]), positions.col(backboneAtoms[0]), positions.col(backboneAtoms[1]), positions.col(backboneAtoms[2]));
            phi_angles[resID - 1] = phi * (180.0 / M_PI);
            phi_indices[resID - 1] = normalizeAngle(phi) / 60.0;
        }
        // Compute psi angles
        if (resID > 0 && resID < numResidues - 1) {
            double psi = computeDihedralAngle(positions.col(backboneAtoms[0]), positions.col(backboneAtoms[1]), positions.col(backboneAtoms[2]), positions.col(backboneAtomsVector[resID+1][0]));
            psi_angles[resID - 1] = psi * (180.0 / M_PI);
            psi_indices[resID - 1] = normalizeAngle(psi)/60.0;
        }
        // Compute plane angles
        if (resID > 0 && resID < numResidues - 2) {
            double plane = computePlaneAngle(positions.col(backboneAtoms[0]), positions.col(backboneAtoms[1]), positions.col(backboneAtoms[3]),
                                             positions.col(backboneAtomsVector[resID+2][0]), positions.col(backboneAtomsVector[resID+2][1]), positions.col(backboneAtomsVector[resID+2][3]));
            plane_angles[resID] = plane * (180.0 / M_PI);
            plane_indices[resID] = normalizeAngle(plane)/30.0;
        }
        // Compute bisect angles
        if (resID > 0 && resID < numResidues - 2) {
            double bisect = computeBisectAngle(positions.col(backboneAtoms[0]), positions.col(backboneAtoms[1]), positions.col(backboneAtoms[3]),
                                               positions.col(backboneAtomsVector[resID+2][0]), positions.col(backboneAtomsVector[resID+2][1]), positions.col(backboneAtomsVector[resID+2][3]));
            bisect_angles[resID] = bisect * (180.0 / M_PI);
            bisect_indices[resID] = normalizeAngle(bisect)/30.0;
        }
    }
}

Eigen::MatrixXd Context::getContacts() const {
    return contacts;
}

Eigen::MatrixXd Context::getBinaryContacts(double cutoff) const {
    return (contacts.array() < cutoff).cast<double>();
}

std::array<Eigen::VectorXi, 4> Context::getAngleIndices() {
    // Get the indices of the angles
    Eigen::VectorXi phi_indices = convertToIndices(60.0, phi_angles);
    Eigen::VectorXi psi_indices = convertToIndices(60.0, psi_angles);
    Eigen::VectorXi plane_indices = convertToIndices(30.0, plane_angles);
    Eigen::VectorXi bisect_indices = convertToIndices(30.0, bisect_angles);
    return {phi_indices, psi_indices, plane_indices, bisect_indices};
}

std::vector<std::array<double, 4>> Context::getSideChainTorsion() {
    return sidechain_torsions;
}

Eigen::VectorXi Context::convertToIndices(double cutoff, Eigen::VectorXd angles) const {
    Eigen::VectorXd indices = angles/cutoff;
    return indices.cast<int>();
}

Eigen::Matrix3Xd Context::getPositions() const {
    return positions;
}

CellList Context::getCellList() const {
    return cellList;
}

Energy Context::getEnergy() const {
    return energy;
}