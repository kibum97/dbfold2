#ifndef HBONDS_H
#define HBONDS_H

#include "globals.h"

#include <string_view>

float HydrogenBonds(
    struct Topology *top,
    struct Context *ctx,
    struct System *sys
);
float FoldHydrogenBonds(
    struct Topology *top,
    struct Context *ctx,
    struct System *sys
);
void InitializeHydrogenBonding(
    struct Topology *top,
    struct Context *ctx,
    struct System *sys,
    struct Simulation *sim
);
int MatchAtomname(std::string_view name);
long int CheckHBond(
    struct Context *ctx,
    struct System *sys,
    struct Topology *top,
    int A, int B);

#endif