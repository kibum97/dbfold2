#ifndef ENERGY_H
#define ENERGY_H

#include "define.h"
#include "globals.h"

float FullAtomEnergy(
    struct System *sys,
    struct Context *ctx
);
void ResetEnergies(
    struct Context *ctx,
    struct Simulation *sim,
    struct System *sys,
    struct Topology *top,
    long int check);
float torsionenergy(
    struct Context *ctx,
    struct Topology *top,
    struct System *sys
);
float sctenergy(
    struct Context *ctx,
    struct System *sys,
    struct Topology *top
);
void aromatic_center(
    struct Context *ctx,
    int res_no, struct vector *V);
float aromatic_plane(
    struct Context *ctx,
    int res_a, int res_b, struct vector *plane_a, struct vector *plane_b);
float aromaticenergy(
    struct Context *ctx,
    struct System *sys,
    struct Topology *top
);

#endif