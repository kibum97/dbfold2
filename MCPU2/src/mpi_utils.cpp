#include "utils/mpi_utils.h"
#include <iostream>

void mpi_initialize(int *argc, char ***argv, int *rank, int *size) {
    MPI_Init(argc, argv);
    MPI_Comm_rank(MPI_COMM_WORLD, rank);
    MPI_Comm_size(MPI_COMM_WORLD, size);
    std::cout << "[MPI Rank " << *rank << "] argc = " << *argc << std::endl;
    for (int i = 0; i < *argc; ++i) {
        std::cout << "[MPI Rank " << *rank << "] argv[" << i << "] = " << (*argv)[i] << std::endl;
    }
}

void mpi_finalize(void) {
    MPI_Finalize();
}