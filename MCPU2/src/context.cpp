#include "context.h"

//TOOD: Need to specify way to define cell size
Context::Context(Topology& topology, Eigen::Matrix3Xd positions)
    : topology(topology), numAtoms(topology.getNumAtoms()), positions(positions), cellList(positions, 5.0, topology.getPeriodicBox()), contacts(0, 0) {
    // Constructor implementation
    contacts = Eigen::MatrixXd::Constant(numAtoms, numAtoms, 10.0);
}

Context::~Context() {
    // Destructor implementation
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

void Context::updateContext(std::vector<size_t> movedAtomIDs, Eigen::Matrix3Xd new_positions) {
    updatePositions(movedAtomIDs, new_positions);
    updateContacts();
}

void Context::updatePositions(std::vector<size_t> movedAtomIDs, Eigen::Matrix3Xd new_positions) {
    this->positions = new_positions;
    this->movedAtomIDs = movedAtomIDs;
}

void Context::updateContacts() {
    std::unordered_set<size_t> movedCellIDs = cellList.updateAtoms(movedAtomIDs, positions);
    for (auto cellID : movedCellIDs) {
        if (cellID >= cellList.cells.size()) exit(1);
        computeContacts(cellID);
    }
}

Eigen::MatrixXd Context::getContacts() const {
    return contacts;
}

Eigen::MatrixXd Context::getBinaryContacts(double cutoff) const {
    return (contacts.array() < cutoff).cast<double>();
}

Eigen::Matrix3Xd Context::getPositions() const {
    return positions;
}

CellList Context::getCellList() const {
    return cellList;
}