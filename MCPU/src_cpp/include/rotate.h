#ifndef ROTATE_H
#define ROTATE_H

#include "globals.h"

#include "atom.h"
#include "vector.h"

void DoRotation(
    struct Context *ctx,
    int a, int b, int c, int d, Float delta_angle, short rotate_natoms,
    short *rotate_atom);

#endif