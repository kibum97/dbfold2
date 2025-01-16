#ifndef MPI_UTIL_H
#define MPI_UTIL_H

#include <mpi.h>

void mpi_initialize(int *argc, char ***argv, int *rank, int *size);
void mpi_finalize(void);

#endif // MPI_UTIL_H