#include "globals.h"

struct Topology     top;
struct System       sys;
struct MCIntegrator integrator;
struct Context      ctx;
struct Simulation   sim;


// --- Rotation and Matrix Variables ---
float rot_mat[3][3];

float rot_mat_00, rot_mat_10, rot_mat_20;
float rot_mat_01, rot_mat_11, rot_mat_21;
float rot_mat_02, rot_mat_12, rot_mat_22;

unsigned char *is_rotated;

// --- Temporary Atom & Coordinate Variables ---
struct atom *temp_atom;
struct atom *temp_prev_atom;
struct int_vector temp_xyz_int; 

// --- Cell / Grid Array Variables ---
struct cell *temp_cell;
struct cell *temp_cell3;
struct cell **temp_cell_array;

// --- Hydrogen Bond & Math Variables ---
float d_memory;
short M, N, O, P;
short *A, *B;

// --- Debugging & Constraints ---
long int distance;
struct alignment constraint_align;