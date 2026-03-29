#include "cell_list.h"
#include <iostream>

CellList::CellList(float cell_size) : inv_cell_size_(1.0f/cell_size) {
    // Constructor implementation
}

CellList::~CellList() {
    // Destructor implementation
}

Eigen::Vector3i CellList::get_cell_index(const Eigen::Vector3f& coords) const {
    // Get the cell ID based on the coordinates
    int cell_x = static_cast<int>(coords(0) * inv_cell_size_);
    int cell_y = static_cast<int>(coords(1) * inv_cell_size_);
    int cell_z = static_cast<int>(coords(2) * inv_cell_size_);
    return Eigen::Vector3i(cell_x, cell_y, cell_z);
}

void CellList::build_cells(const Eigen::Matrix3Xf& positions) {
    // Build cells based on the positions
    // This is a placeholder implementation
    int num_atoms = positions.cols();
    for (int i = 0; i < num_atoms; ++i) {
        Eigen::Vector3i cell_id = get_cell_index(positions.col(i));
        cells[cell_id].atom_indices.push_back(i);
    }
}

void CellList::update_moved_atoms(std::vector<int>& moved_atoms, Eigen::Matrix3Xf& old_positions, Eigen::Matrix3Xf& new_positions) {
    // Update the list of moved atoms
    // This is a placeholder implementation
    for (int atom_id : moved_atoms) {
        Eigen::Vector3i old_cell = get_cell_index(old_positions.col(atom_id));
        Eigen::Vector3i new_cell = get_cell_index(new_positions.col(atom_id));

        if (old_cell != new_cell) {
            // Remove from old cell
            auto& old_vec = cells[old_cell].atom_indices;
            old_vec.erase(std::remove(old_vec.begin(), old_vec.end(), atom_id), old_vec.end());
            // Add to new cell
            cells[new_cell].atom_indices.push_back(atom_id);
        }
    }
}

std::vector<int> CellList::get_neighbor_atom_indices(const Eigen::Vector3f& pos) const {
    std::vector<int> result;
    Eigen::Vector3i center = get_cell_index(pos);

    for (int dx = -1; dx <= 1; ++dx) {
        for (int dy = -1; dy <= 1; ++dy) {
            for (int dz = -1; dz <= 1; ++dz) {
                Eigen::Vector3i neighbor = center + Eigen::Vector3i(dx, dy, dz);
                auto it = cells.find(neighbor);
                if (it != cells.end()) {
                    result.insert(result.end(), it->second.atom_indices.begin(), it->second.atom_indices.end());
                }
            }
        }
    }
    return result;
}