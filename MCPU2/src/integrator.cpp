#include "integrator.h"
#include "utils/random.h"

Integrator::Integrator() {
    // Constructor
}

Integrator::~Integrator() {
    // Destructor
}

Eigen::Matrix3d Integrator::createRotationMatrix(double angle, Eigen::Vector3d axis) {
    // Create a rotation matrix
    Eigen::Matrix3d rotationMatrix;
    rotationMatrix = Eigen::AngleAxisd(angle, axis);
    return rotationMatrix;
}

Eigen::Quaterniond Integrator::createQuaternion(double angle, Eigen::Vector3d axis) {
    return Eigen::Quaterniond(Eigen::AngleAxisd(angle, axis));
}

Eigen::Matrix3Xd Integrator::partialRotation(Eigen::Matrix3Xd& positions, const std::vector<size_t>& atomIDs, double angle, Eigen::Vector3d axis) {
    // Set origin
    Eigen::Vector3d origin = positions.col(atomIDs[1]);
    // Translate the molecule so that atom 2 is at the origin
    for (int i = 0; i < positions.cols(); ++i) {
        positions.col(i) -= origin;
    }
    // Create a rotation matrix
    Eigen::Matrix3d rotationMatrix = createRotationMatrix(angle, axis);
    // Eigen::Quaterniond quaternion = createQuaternion(angle, axis);
    for (auto atomid = atomIDs.begin() + 2; atomid != atomIDs.end(); ++atomid) {
        size_t atomID = *atomid;
        positions.col(atomID) = rotationMatrix * positions.col(atomID);
    }
    // Translate the molecule back to its original position
    for (int i = 0; i < positions.cols(); ++i) {
        positions.col(i) += origin;
    }
    return positions;
}

std::vector<size_t> Integrator::selectAtomIDs(size_t atomID1, size_t atomID2, const Topology& topology) {
    // Select the atoms that moves
    std::vector<size_t> atomIDs;
    // For the residue that the axis belong to
    // we have to select subset of atoms based on the axis
    std::string atomName2 = topology.atoms[atomID2].atomName;
    for (const auto& atom : topology.residues[topology.atoms[atomID2].resID].atoms) {
        if (atomName2 == "N") {
            atomIDs.push_back(atom.atomID);
        } else if (atomName2 == "CA") {
            if (atom.atomName != "N" && atom.atomName != "H") {
                atomIDs.push_back(atom.atomID);
            }
        } else if (atomName2 == "C") {
            if (atom.atomName == "C" | atom.atomName == "O") {
                atomIDs.push_back(atom.atomID);
            }
        } else {
            std::cerr << "Error: Invalid atom for rotation" << std::endl;
            exit(EXIT_FAILURE);
        }
    }
    // For the remaining residues
    // add all the atoms that belong to them
    for (size_t residueID = topology.atoms[atomID2].resID + 1; residueID <= topology.residues.size(); ++residueID) {
        for (const auto& atom : topology.residues[residueID].atoms) {
            atomIDs.push_back(atom.atomID);
        }
    }
    return atomIDs;
}

Eigen::Matrix3Xd Integrator::sidechainRotation(size_t resID, Eigen::Matrix3Xd& positions, std::unordered_map<size_t, std::vector<RotamerData>> rotamerIDMap, std::unordered_map<size_t, std::vector<TorsionIDData>> rotatingSCAtomsMap, std::vector<std::array<double, 4>> sidechain_torsions) {
    // Get rotation angles from rotamer data
    std::cout << "Performing sc rotation at residue " << resID << std::endl;
    std::vector<size_t> rotamerIDs = {};
    std::vector<double> weights = {};
    size_t index = 0;
    for (const auto& rotamerData : rotamerIDMap[resID]) {
        rotamerIDs.push_back(index);
        weights.push_back(rotamerData.probability);
        ++index;
    }
    if (rotamerIDs.size() > 0) {
        RotamerData selected_rotamer = rotamerIDMap[resID][selectElementWithWeights(rotamerIDs, weights)];
        std::cout << "Selected rotamer confirmed" << std::endl;
        // Perform sidechain rotation
        for (const auto& torsionIDData : rotatingSCAtomsMap[resID]) {
            int torsion_id = torsionIDData.torsion_id;
            std::vector<size_t> atomIDs;
            atomIDs.push_back(torsionIDData.torsion_atomIDs[1]);
            atomIDs.push_back(torsionIDData.torsion_atomIDs[2]);
            for (size_t atomID : torsionIDData.rotating_atomIDs) {
                atomIDs.push_back(atomID);
            }
            std::cout << "Atom IDs confirmed" << std::endl;
            Eigen::Vector3d axis = positions.col(atomIDs[1]) - positions.col(atomIDs[0]);
            axis.normalize();
            std::cout << "Axis computed" << std::endl;
            double angle = generateGaussian(selected_rotamer.chiMeans[torsion_id],selected_rotamer.chiStdDevs[torsion_id]) - sidechain_torsions[resID][torsion_id];
            std::cout << "Angle computed" << std::endl;
            positions = partialRotation(positions, atomIDs, angle, axis);
        }
    }
    return positions;
}

/*
void Integrator::concertedRotation(Eigen::Matrix3Xd& positions, const std::vector<size_t>& atomIDs, double angle    ) {
    // Perform concerted rotation
    std::vector<size_t> movedAtomIDs;
    // Get the atom's residue and its atoms
    size_t residueID = topology.getResidueID(atomID);
    std::vector<size_t> residueAtomIDs;
    for (auto i = 0; i < 3; ++i) {
        for (const auto& atom : topology.residues[residueID + i].atoms) {
            movedAtomIDs.push_back(atom.atomID);
        }
    }
    
    
    topology.residues[residueID].atoms

}
*/