#include "globals.h"

/* MPI and System */
int nprocs = 0;
int myrank = 0;
int natoms = 0;
int nresidues = 0;
FILE *STATUS = NULL;
MPI_Comm mpi_world_comm;

/* Pointers to Data Structures */
struct atom *native = NULL;
struct atom *prev_native = NULL;
struct atom *orig_native = NULL;
struct atom *native_Emin = NULL;
struct atom *native_RMSDmin = NULL;
struct residue *native_residue = NULL;
struct contact_data **data = NULL;
struct contact_data **struct_data = NULL;
struct cell ***the_matrix = NULL;
struct pair *ab = NULL;
struct pair *cd = NULL;

/* Simulation State */
struct monte_carlo mc;
Float rot_mat[3][3];
short nclashes = 0;
short ncontacts = 0;
long int mcstep = 0;
long int seed = 0;

/* Energy Variables */
Float E = 0, E_pot = 0, dE = 0, dE_pot = 0, prev_E = 0, prev_E_pot = 0;
Float E_tor = 0, E_sct = 0, E_aro = 0, E_Rg = 0;
Float dE_tor = 0, dE_sct = 0, dE_aro = 0, dE_Rg = 0;
float Emin = 0;

/* Constraint Globals */
int Max_N_constraints = 1000;
int N_constraints = 0;
int **constraint_array = NULL;
float *distances = NULL;
float *constraint_weights = NULL;
float k_constraint = 0.0;
float E_constraint = 0.0;
float prev_E_constraint = 0.0;
float dE_constraint = 0.0;
float mean_constraint_distance = 0.0;

/* Replica and MPI Buffers */
float *Enode = NULL, *Tnode = NULL, *buf_in = NULL, *buf_out = NULL;
int *Nnode = NULL, *Cnode = NULL;