#define FORCEFIELD 'MCPU'

#include <stdio.h>
#include <stdlib.h>
#include "mpi_util.h"
#include "read_config.h"
#include "topology.h"
#include "logger.h"
#include "load_pdb.h"
#include "context.h"

int main(int argc, char *argv[]) {
    int rank, size;

    // Initialize MPI
    mpi_initialize(&argc, &argv, &rank, &size);

    // Initialize Logger
    if (rank == 0) {
        Logger::initialize("simulation.log");
        Logger::log(INFO, "Simulation started");
    }

    // Check if configuration file is provided
    if (argc < 2) {
        if (rank == 0) {
            Logger::logf(ERROR, "Usage: %s <config_file>", argv[0]);
        }
        mpi_finalize();
        return EXIT_FAILURE;
    }

    // Read configuration file
    if (rank == 0) {
        SetProgramOptions(argv[1]);
        LogConfigurationDetails();
        
        // Initialize Topology and Forcefield
        std::cout << "Parsing PDB file..." << std::endl;
        auto [topology, positions] = parsePDB(native_file);
        std::cout << "Parsed PDB file successfully." << std::endl;

        // Initialize Context
        std::cout << "Initializing context..." << std::endl;
        Context context(topology, positions);
        std::cout << "Context initialized successfully." << std::endl;

        CellList cellList = context.getCellList();
        for (int i = 0; i < cellList.cells.size(); ++i) {
            std::cout << "Atom IDs in cell " << i << ": ";
            for (const auto& atomID : cellList.cells[i].atomIDs) {
                std::cout << atomID << " ";
            }
            std::cout << std::endl;
        }
    }

    // Initialize matrix

    // Finalize MPI
    mpi_finalize();

    // Finalize Logger
    if (rank == 0) {
        Logger::log(INFO, "Simulation ended");
        Logger::finalize();
    }

    return EXIT_SUCCESS;
}