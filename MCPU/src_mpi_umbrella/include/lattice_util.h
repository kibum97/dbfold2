#ifndef LATTICE_UTIL_H
#define LATTICE_UTIL_H

#include "globals.h" /* for struct atom */

struct atom;

/* MATRIX_SIZE is set to be 20 in init.h (it's a u. char, not a macro) */

void          FindLatticeCoordinates(struct atom *);
void          UpdateLattice(short, short *);
void          InitializeMatrix();
unsigned char PerBound(signed char);
void          CopyLatticeCoordinates(struct atom, struct atom *);

#endif