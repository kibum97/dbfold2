#include "globals.h"

#include "globals.h"

/* Define and initialize variables here ONCE */
int Max_N_constraints = 1000;
int **constraint_array = NULL;
float *distances = NULL;
float *constraint_weights = NULL;
int *disulfide_pairs_attempt = NULL;
float k_constraint = 0.0;
float mean_constraint_distance = 0.0;
float CLUSTER_MOVE = 0.33;
float USE_CLUSTER = 0.0;
int MAX_CLUSTERSTEP = 0;
int move_cycle = 0;
char amino_letters[20] = "ARNDCQEGHILKMFPSTWYV";

/* MPI and Logging Globals */
int nprocs;
int myrank;
MPI_Comm mpi_world_comm;
FILE *STATUS = NULL;

/* Matrix Pointer */
struct cell ***the_matrix = NULL;