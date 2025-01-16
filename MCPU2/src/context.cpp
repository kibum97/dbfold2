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
            std::cout << "Atom ID: " << atomID1 << std::endl;
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

void Context::updatePositions(std::map<size_t, Eigen::Vector3d> new_positions) {
    for (auto const& [atom_id, new_position] : new_positions) {
        if (atom_id >= positions.cols()) exit(1);
        positions.col(atom_id) = new_position;
    }
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

std::vector<size_t> Context::setMovedAtomIDs() {
    return movedAtomIDs;
}