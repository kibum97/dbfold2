#ifndef DEBUG_H
#define DEBUG_H

#if DEBUG
#include "globals.h"

void Debug(
    struct Context *ctx,
    struct Topology *top,
    struct Simulation *sim,
    struct System *sys
);
void DebugContacts(
    struct Simulation *sim,
    struct Topology *top,
    struct Context *ctx,
    struct System *sys
);
void CheckForDebugContacts(
    struct Simulation *sim,
    struct Topology *top,
    struct Context *ctx,
    struct System *sys,
    int a, int b
);
#endif

#endif