#ifndef CELL_LIST_H
#define CELL_LIST_H

#include <Eigen/Dense>
#include <array>
#include <map>
#include <vector>

struct Vector3iHash {
    std::size_t operator()(const Eigen::Vector3i &key) const {
        std::size_t h1 = std::hash<int>{}(key[0]);
        std::size_t h2 = std::hash<int>{}(key[1]);
        std::size_t h3 = std::hash<int>{}(key[2]);
        return h1 ^ (h2 << 1) ^ (h3 << 2);
    }
};

struct Cell {
    std::vector<int> atom_indices;
};

class CellList {
   public:
    CellList(float cell_size);
    ~CellList();

    std::unordered_map<Eigen::Vector3i, Cell, Vector3iHash> cells;
    std::map<int, int>                                      atom_to_cell_map;

    Eigen::Vector3i  get_cell_index(const Eigen::Vector3f &coords) const;
    void             build_cells(const Eigen::Matrix3Xf &positions);
    void             update_moved_atoms(std::vector<int> &moved_atoms, Eigen::Matrix3Xf &positions,
                                        Eigen::Matrix3Xf &new_positions);
    std::vector<int> get_neighbor_atom_indices(const Eigen::Vector3f &pos) const;

   private:
    float inv_cell_size_;
};

#endif