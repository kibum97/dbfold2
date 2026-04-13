#include "CellGrid.hpp"
#include <Eigen/Dense>
#include <algorithm>

// Note: Assuming coords is Matrix4Xf where rows 0,1,2 are X,Y,Z and row 3 is W (charge/mass/pad)
void BuildDynamicGrid(CellGrid& grid, const Eigen::Matrix4Xf& coords, double cell_width) {
    grid.num_atoms = coords.cols(); // Number of columns = number of atoms
    grid.inverse_size = 1.0 / cell_width;

    // 1. Vectorized Bounding Box
    // topRows<3>() strictly looks at X, Y, Z, ignoring the 4th row.
    // rowwise().minCoeff() executes heavily optimized SIMD reductions.
    Eigen::Vector3f min_c = coords.topRows<3>().rowwise().minCoeff();
    Eigen::Vector3f max_c = coords.topRows<3>().rowwise().maxCoeff();

    // Offset origin slightly
    grid.origin_x = min_c.x() - 0.001;
    grid.origin_y = min_c.y() - 0.001;
    grid.origin_z = min_c.z() - 0.001;

    // 2. Calculate dynamic grid dimensions
    grid.size_x = (int)((max_c.x() - grid.origin_x) * grid.inverse_size) + 1;
    grid.size_y = (int)((max_c.y() - grid.origin_y) * grid.inverse_size) + 1;
    grid.size_z = (int)((max_c.z() - grid.origin_z) * grid.inverse_size) + 1;

    int total_cells = grid.size_x * grid.size_y * grid.size_z;

    // 3. Grow vectors only if necessary
    if (grid.head.size() < total_cells) {
        grid.head.resize(total_cells);
    }
    if (grid.next.size() < grid.num_atoms) {
        grid.next.resize(grid.num_atoms);
        grid.atom_to_cell.resize(grid.num_atoms);
    }

    // 4. Zero out active cells
    std::fill(grid.head.begin(), grid.head.begin() + total_cells, -1);

    // 5. Rebuild the flat Head/Next linked cell list
    // A single pass using Eigen's fast coeff-access
    for (int i = 0; i < grid.num_atoms; ++i) {
        // coords(row, col) is fast and contiguous in column-major memory
        int gx = (int)((coords(0, i) - grid.origin_x) * grid.inverse_size);
        int gy = (int)((coords(1, i) - grid.origin_y) * grid.inverse_size);
        int gz = (int)((coords(2, i) - grid.origin_z) * grid.inverse_size);
        
        int cell_id = (gx * grid.size_y * grid.size_z) + (gy * grid.size_z) + gz;
        
        grid.next[i] = grid.head[cell_id];
        grid.head[cell_id] = i;
        grid.atom_to_cell[i] = cell_id;
    }
}