#ifndef ATOM_H
#define ATOM_H
#include "globals.h" /* for global variables */

struct low_resol {
    struct vector xyz;
};

void CopyAtom(struct atom, struct atom *);

#endif
