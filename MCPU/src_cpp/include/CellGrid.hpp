#include <vector>
#include <cmath>

struct CellGrid {
    int size;             
    double inverse_size;  
    int half_size;        
    int num_atoms;
    
    // Pure flat arrays. Pre-allocated ONCE at startup.
    std::vector<int> head;        // Size: total number of cells. 
    std::vector<int> next;        // Size: total number of atoms.
    std::vector<int> atom_to_cell; // Size: total number of atoms. Maps atom -> current cell
    
    // Cell neighbor map (flat 2D representation)
    // Size: total_cells * 27. Stores 1D indices of neighboring cells.
    std::vector<int> cell_neighbors; 
};