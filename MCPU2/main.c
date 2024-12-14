#include <stdio.h>
#include <stdlib.h>
#include "mpi_util.h"

int main(int argc, char *argv[]) {
    int rank, size;
    int data, received_data;

    // Initialize MPI
    mpi_initialize(&argc, &argv, &rank, &size);

    // Example: Broadcast data from rank 0 to all other ranks
    if (rank == 0) {
        data = 100; // Data to be broadcasted
    }
    MPI_Bcast(&data, 1, MPI_INT, 0, MPI_COMM_WORLD);
    printf("Rank %d received data: %d\n", rank, data);

    // Example: Send and receive data between ranks
    if (rank == 0) {
        data = 200; // Data to be sent
        MPI_Send(&data, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
        printf("Rank %d sent data: %d\n", rank, data);
    } else if (rank == 1) {
        MPI_Recv(&received_data, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("Rank %d received data: %d\n", rank, received_data);
    }

    // Example: Reduce data from all ranks to rank 0
    data = rank; // Each rank contributes its rank number
    int sum;
    MPI_Reduce(&data, &sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        printf("Sum of ranks: %d\n", sum);
    }

    // Finalize MPI
    mpi_finalize();

    return 0;
}