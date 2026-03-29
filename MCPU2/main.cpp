#define FORCEFIELD 'MCPU'
#include <stdio.h>
#include <stdlib.h>
#include "utils/mpi_utils.h"
#include <iomanip>
#include "logger.h"
#include "read_config.h"
#include "pdb_parser.h"

#include "forces/ForcefieldMCPU.h"
#include "system_mcpu.h"
#include "context_mcpu.h"

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstddef>

/*
TODO
- implement hydrogen bond energy and aromatic energy and torsional energy
- implement the move acceptance criteria
- implement the MC steps
- implement the replica exchange
- Make reading files part in the same style
*/


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
        std::cout << "Usage: " << argv[0] << " <config_file>" << std::endl;
        if (rank == 0) {
            Logger::logf(ERROR, "Usage: %s <config_file>", argv[0]);
        }
        mpi_finalize();
        return EXIT_FAILURE;
    }

    if (rank == 0) {
        // Read configuration file
        SetProgramOptions(argv[1]);
        LogConfigurationDetails();

        // Initialize Topology
        // At this stage topology will be simply read from a PDB file
        // without considering the forcefield
        std::cout << "Parsing PDB file..." << std::endl;
        PDBFile pdbfile(native_file);
        pdbfile.protein.validate_structure();
        pdbfile.protein.print_structure();
        std::cout << "Parsed PDB file successfully. Make sure that you are using the PDB without hydrogen atoms" << std::endl;

        // Initialize Context and System
        // This is where the forcefield dependent parameters will be set
        std::cout << "Initializing System..." << std::endl;
        ForcefieldMCPU forcefield;
        SystemMCPU system(forcefield, pdbfile.protein);
        std::cout << "System initialized successfully" << std::endl;
        // Initialize Context
        std::cout << "Initializing Context..." << std::endl;
        ContextMCPU context(system, pdbfile.positions);
        std::cout << "Context initialized successfully" << std::endl;
        //context.compute_energy();
        std::cout << "Energy computed successfully" << std::endl;
        // Print energy values
        std::cout << "Mu Energy: " << context.energy.mu_energy << std::endl;
        std::cout << "Backbone Torsion Energy: " << context.energy.bb_torsion_energy << std::endl;
        std::cout << "Sidechain Torsion Energy: " << context.energy.sc_torsion_energy << std::endl;
        std::cout << "Aromatic Energy: " << context.energy.aromatic_energy << std::endl;
        std::cout << "Hydrogen Bond Energy: " << context.energy.hbond_energy << std::endl;
        std::cout << "Total Energy: " << context.energy.total_energy << std::endl;

        // DEBUG
        for (const auto& chain: pdbfile.protein.chains) {
            for (const auto& residue : chain->residues) {
                int chain_id = chain->id;
                int residue_id = residue->id;
                int global_index = system.get_global_residue_index(chain_id, residue_id);
            }
        }
        std::cout << "Done" << std::endl;
    }

    // Finalize Logger
    if (rank == 0) {
        std::cout << "Finalizing Logger..." << std::endl;
        Logger::log(INFO, "Simulation ended");
        Logger::finalize();
    }

    // Finalize MPI
    mpi_finalize();
    std::cout << "Exiting safely." << std::endl;
    return 0;
}