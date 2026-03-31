#include <mpi.h>
#include <stdio.h>

#include "constraint.h"
#include "fold.h"
#include "globals.h"
#include "init.h"
#include "lattice_util.h"
#include "loop.h"
#include "setup.h"

/*
// 1. Build the static structural definition
Topology topology = PDBParser::Load("protein.pdb");

// 2. Define the physics and force weights
System system;
system.addForce(new HydrogenBondForce(weight_hbond, SEQ_DEP_HB));
system.addForce(new ClashForce(weight_clash, NON_NATIVE_REPULSION));

// 3. Define how the system moves (The MC variables)
MonteCarloIntegrator integrator(MC_TEMP, STEP_SIZE, CLUSTER_MOVE);

// 4. Create the living state
Context context(system, integrator, topology);

// 5. Wrap it all in a Simulation and attach I/O
Simulation simulation(topology, system, integrator, context);
simulation.addReporter(new PDBReporter("output.pdb", MC_PDB_PRINT_STEPS));
simulation.addReporter(new LogReporter("status.log", MC_PRINT_STEPS));

// Run the engine
simulation.step(MC_STEPS);
*/

int main(int argc, char *argv[]) {
    /* Input safety check */
    if (argc != 2) {
        printf("ERROR! Usage: ./fold_potential config_file\nargc : %d\n", argc);
        exit(1);
    }

    /* Initialize simulation defaults */
    BootSimulationDefaults();

    /* MPI initialization */
    if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
        fprintf(stderr, "Error initializing MPI.\n");
        return 1;
    } else if (MPI_Comm_size(MPI_COMM_WORLD, &sim.nprocs) != MPI_SUCCESS) {
        fprintf(stderr, "Error getting MPI communicator size.\n");
        return 1;
    } else if (MPI_Comm_rank(MPI_COMM_WORLD, &sim.myrank) != MPI_SUCCESS) {
        fprintf(stderr, "Error getting MPI rank.\n");
        return 1;
    }

    /* current_replica tracks a set of coordinates across replica space */
    sim.current_replica = sim.myrank;

    /* initialize data */
    SetProgramOptions(&sim, &sys, &integrator, &ctx, &top, argc, argv);

    sim.seed = time(NULL);
    sim.seed += (int)(integrator.MC_TEMP * 1000) + getpid() + sim.myrank;
    // srand48(seed);
    set_threefry_array((unsigned long int)&sim.seed);
    // fprintf(STATUS, "nprocs = %d, myrank = %d\n", nprocs, myrank);
    fprintf(sim.STATUS, "---GENERAL---\n");
    fprintf(sim.STATUS, "  seed:\t\t%ld\n\n", sim.seed);
    fprintf(sim.STATUS, "  VERSION: flat-bottom harmonic disulfide potential for all distances\n");
    fflush(sim.STATUS);

    sys.weight_potential = POTNTL_WEIGHT;
    sys.weight_hbond     = HBOND_WEIGHT;

    integrator.STEP_SIZE = 2 * deg2rad;

    InitializeMatrix(&ctx, &integrator); /*AB: Initializes a matrix whose function I'm not totally sure of at the
                           moment*/
    InitializeProtein(&ctx, &top, &sim, &sys, &integrator);
    if (strcmp(sim.rmsd_constraint_file, "None") != 0) {
        Init_frag_rmsd_constraints(&sim);
    }
    if (strcmp(sim.movable_region_file, "None") != 0) {  // KP added for movable region constraint
        Read_movable_region(&sim, &top);
    } else {
        top.n_movable_residues = top.nresidues;
        for (int i = 0; i < top.nresidues; i++) {
            top.movable_residue_map[i] = i;
        }
    }
    for (int i = 0; i < top.nresidues; i++)
        top.total_ntorsions += ctx.native_residue[i].ntorsions;

    if (integrator.CLUSTER_MOVE && integrator.USE_CLUSTER) {
        fprintf(sim.STATUS,
                "WARNING! KNOWLEDGE MOVES ARE IMPLEMENTED! THESE DO NOT SATISFY DETAILED BALANCE!  "
                "\n ");
    }
    Fold(&ctx, &top, &sim, &sys, &integrator);

    /* vz: We ought to free memory we allocated, but I guess the program
     * will exit soon, and during program run, there aren't memory
     * leaks. */

    MPI_Finalize();
    return 0;
}