#ifndef UPDATE_H
#define UPDATE_H

#include "globals.h"

void Restore(
    struct Context *ctx,
    struct MCIntegrator *integrator,
    struct System  *sys
);
void Update(
    struct Context *ctx,
    struct System  *sys,
    struct Simulation *sim
);

#endif