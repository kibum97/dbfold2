#ifndef FOLD_H
#define FOLD_H

#include <time.h>

#include "globals.h"

void SetupMatrixStuff(
    struct MCIntegrator *integrator,
    struct Topology *top,
    struct Context *ctx
);
void Fold(
    struct Context *ctx,
    struct Topology *top,
    struct Simulation *sim,
    struct System *sys,
    struct MCIntegrator *integrator
);

#endif