#include "cell_list.h"
#include <algorithm>
#include <iostream>

void Cell::addAtomID(size_t atomID) {
    atomIDs.push_back(atomID);
}

void Cell::removeAtomID(size_t atomID) {
    atomIDs.erase(std::remove(atomIDs.begin(), atomIDs.end(), atomID), atomIDs.end());
}

CellList::CellList(Eigen::Matrix3Xd positions, double cellSize, const Eigen::Matrix3Xd& periodicBox) : cellSize(cellSize) {
    Eigen::Vector3d boxSize = periodicBox.diagonal();
    numCells = {static_cast<int>(boxSize[0] / cellSize),
                static_cast<int>(boxSize[1] / cellSize),
                static_cast<int>(boxSize[2] / cellSize)};
    cells.resize(numCells[0] * numCells[1] * numCells[2]);
    for (int i = 0; i < numCells[0]; ++i) {
        for (int j = 0; j < numCells[1]; ++j) {
            for (int k = 0; k < numCells[2]; ++k) {
                std::array<int, 3> index = {i, j, k};
                size_t flatIndex = flattenIndex(index);
                cells[flatIndex].cellID = flatIndex;
                cells[flatIndex].neighborIds = mapNeighborIds(index);
                cells[flatIndex].atomIDs.reserve(static_cast<int>(cellSize * cellSize * cellSize));
            }
        }
    }
    for (size_t atomID = 0; atomID < positions.cols(); ++atomID) {
        std::array<int, 3> cellIndex = getCellIndex(positions.col(atomID));
        size_t cellId = flattenIndex(cellIndex);
        cells[cellId].addAtomID(atomID);
    }
    mapCellIds();
}

std::unordered_set<size_t> CellList::updateAtoms(std::vector<size_t>& movedAtomIDs, Eigen::Matrix3Xd positions) {
    std::unordered_set<size_t> movedCellIDs;
    for (auto atomID : movedAtomIDs) {
        int old_cell_id = cellIdMap[atomID];
        std::array<int, 3> cell_index = getCellIndex(positions.col(atomID));
        int cell_id = flattenIndex(cell_index);
        if (old_cell_id != cell_id) {
            cells[old_cell_id].removeAtomID(atomID);
            cells[cell_id].addAtomID(atomID);
        }
        movedCellIDs.insert(old_cell_id);
        movedCellIDs.insert(cell_id);
    }
    return movedCellIDs;
}

std::array<int, 3> CellList::getCellIndex(const Eigen::Vector3d& position) const {
    std::array<int, 3> index = {
        static_cast<int>(position[0] / cellSize),
        static_cast<int>(position[1] / cellSize),
        static_cast<int>(position[2] / cellSize)
    };
    for (int i = 0; i < 3; ++i) {
        index[i] = (index[i] % numCells[i] + numCells[i]) % numCells[i];
    }
    return index;
}

size_t CellList::flattenIndex(const std::array<int, 3>& index) const {
    return static_cast<size_t>(index[0] + numCells[0] * (index[1] + numCells[1] * index[2]));
}

std::vector<size_t> CellList::mapNeighborIds(const std::array<int, 3>& index) const {
    std::vector<size_t> neighborIds;
    std::array<int, 3> neighborIndex;
    neighborIds.reserve(27);
    for (int i = -1; i <= 1; ++i) {
        for (int j = -1; j <= 1; ++j) {
            for (int k = -1; k <= 1; ++k) {
                neighborIndex = {index[0] + i, index[1] + j, index[2] + k};
                for (int l = 0; l < 3; ++l) {
                    if (neighborIndex[l] < 0) {
                        neighborIndex[l] += numCells[l];
                    } else if (neighborIndex[l] >= numCells[l]) {
                        neighborIndex[l] -= numCells[l];
                    }
                }
                neighborIds.push_back(flattenIndex(neighborIndex));
            }
        }
    }
    return neighborIds;
}

void CellList::mapCellIds() {
    for (int i = 0; i < cells.size(); ++i) {
        for (auto atomID : cells[i].atomIDs) {
            cellIdMap[atomID] = i;
        }
    }
}


/*
#include "cell_list.h"

void InitializeCellList(Matrix3D& cell_list, int MATRIX_SIZE) {
    // Resize the 3D matrix
    cell_list.resize(MATRIX_SIZE, std::vector<std::vector<Cell>>(MATRIX_SIZE, std::vector<Cell>(MATRIX_SIZE)));

    for (int i = 0; i < MATRIX_SIZE; ++i) {
        for (int j = 0; j < MATRIX_SIZE; ++j) {
            for (int k = 0; k < MATRIX_SIZE; ++k) {
                cell_list[i][j][k].natoms = 0;
                cell_list[i][j][k].X = i;
                cell_list[i][j][k].Y = j;
                cell_list[i][j][k].Z = k;

                for (int a = -1; a <= 1; ++a) {
                    for (int b = -1; b <= 1; ++b) {
                        for (int c = -1; c <= 1; ++c) {
                            int ni = (i + a + MATRIX_SIZE) % MATRIX_SIZE;
                            int nj = (j + b + MATRIX_SIZE) % MATRIX_SIZE;
                            int nk = (k + c + MATRIX_SIZE) % MATRIX_SIZE;
                            cell_list[i][j][k].neighbors.push_back(&cell_list[ni][nj][nk]);
                        }
                    }
                }
            }
        }
    }
}

size_t PerBound(int X, int MATRIX_SIZE) {
    return (X + MATRIX_SIZE) % MATRIX_SIZE;  // Handles underflow and overflow
}
*/