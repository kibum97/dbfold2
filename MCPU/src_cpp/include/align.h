#ifndef ALIGN_H
#define ALIGN_H

#include "globals.h"

void ReadAlignment(
    struct Topology *top,
    struct System *sys
);
void ReadHelicityData(
    struct System *sys,
    struct Topology *top,
    struct Context *ctx,
    struct Simulation *sim
);
void SetupAlignmentStructure(
    struct Context *ctx,
    struct Topology *top,
    struct System *sys,
    struct MCIntegrator *integrator,
    struct Simulation *sim
);
void SetupAlignmentPotential(
    struct Topology *top,
    struct Context *ctx,
    struct System *sys,
    struct Simulation *sim
);

#endif