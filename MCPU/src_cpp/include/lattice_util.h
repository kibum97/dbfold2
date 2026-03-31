#ifndef LATTICE_UTIL_H
#define LATTICE_UTIL_H

#include "globals.h" /* for struct atom */

struct atom;

/* MATRIX_SIZE is set to be 20 in init.h (it's a u. char, not a macro) */

void FindLatticeCoordinates(
    struct Context *ctx, struct MCIntegrator *integrator,
    struct atom *ATOM);
void UpdateLattice(
    struct Context *ctx, struct MCIntegrator *integrator,
    short rotate_natoms, short *rotate_atom);
void InitializeMatrix(
    struct Context *ctx, struct MCIntegrator *integrator
);
unsigned char PerBound(
    struct MCIntegrator *integrator,
    signed char X);
void CopyLatticeCoordinates(struct atom A, struct atom *B);

#endif