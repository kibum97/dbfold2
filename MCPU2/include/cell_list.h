#ifndef CELL_LIST_H
#define CELL_LIST_H

#include <vector>
#include <array>
#include <Eigen/Dense>
#include <unordered_map>
#include <map>
#include <unordered_set>

class Cell {
public:
    size_t cellID = 0;
    std::vector<size_t> atomIDs = {};
    std::vector<size_t> neighborIds = {};
    // Functions handling atomIDs
    void addAtomID(size_t atomID);
    void removeAtomID(size_t atomID);
};

class CellList {
public:
    CellList(Eigen::Matrix3Xd positions, double cellSize, const Eigen::Matrix3Xd& periodicBox);

    std::vector<Cell> cells;
    std::array<int, 3> numCells;
    std::map<int, int> cellIdMap; // Map atom id to cell id
    double cellSize;

    std::unordered_set<size_t> updateAtoms(std::vector<size_t>& movedAtomIDs, Eigen::Matrix3Xd positions);
    std::array<int, 3> getCellIndex(const Eigen::Vector3d& position) const;
    size_t flattenIndex(const std::array<int, 3>& index) const;
    std::vector<size_t> mapNeighborIds(const std::array<int, 3>& index) const;
    void mapCellIds();
};

#endif // CELL_LIST_H






/*
#ifndef CELL_LIST_H
#define CELL_LIST_H

#include <Eigen/Dense>
#include <vector>
#include "constants.h"

struct Cell {
    int natoms;
    int X, Y, Z;
    std::vector<Cell*> neighbors; // Pointer to neighboring cells
};

using Matrix3D = std::vector<std::vector<std::vector<Cell>>>;

void InitializeCellList(Matrix3D& cell_list, int MATRIX_SIZE);

#endif // CELL_LIST_H
*/