#include "mpi_util.h"

void mpi_initialize(int *argc, char ***argv, int *rank, int *size) {
    MPI_Init(argc, argv);
    MPI_Comm_rank(MPI_COMM_WORLD, rank);
    MPI_Comm_size(MPI_COMM_WORLD, size);
}

void mpi_finalize(void) {
    MPI_Finalize();
}