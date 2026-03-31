#ifndef CONSTRAINT_H
#define CONSTRAINT_H

#include "getrms.h"
#include "globals.h"

extern struct alignment constraint_align;

void Init_frag_rmsd_constraints(
    struct Simulation *sim
);
float Compute_constraint_energy(
    struct Context *ctx, struct System *sys,
    const struct residue *residues, const struct atom *atoms);
float Compute_frag_rmsd_constraint_energy(
    struct System *sys,
    struct backbone  struct1[MAXSEQUENCE],
    struct backbone  struct2[MAXSEQUENCE],
    struct alignment algn);
int Argmin(
    struct System *sys,
    float distances[]);
void Read_constraints(
    struct System *sys,
    struct Simulation *sim,
    struct Context *ctx
);

#endif