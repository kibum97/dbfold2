#include "globals.h"

void BootSimulationDefaults() {
    /* 1. System Physics Defaults */
    sys.contact_calpha_cutoff = 7.0;
    sys.k_bias                = 0.0;
    sys.USE_GO_POTENTIAL      = 0;
    sys.MAX_TYPES             = 20;  // Assuming this is your standard

    /* 2. Integrator Defaults */
    integrator.USE_CLUSTER = 0;
    integrator.MC_TEMP     = 1.0;

    /* 3. Context Defaults (The starting state) */
    ctx.rms_RMSDmin    = 100.0;
    ctx.rms_Emin       = 100.0;
    ctx.mcstep_RMSDmin = 0;
    ctx.mcstep_Emin    = 0;
    ctx.native_rms     = 0.0;

    /* 4. Simulation / I/O Defaults */
    sim.fpdb   = NULL;
    sim.STATUS = NULL;
    sim.DATA   = NULL;

    // WARNING: In C, you cannot use `=` to assign strings after creation!
    // You must use strcpy or snprintf.
    strcpy(sim.path_dir, "./");
    strcpy(sim.native_file, "");
}