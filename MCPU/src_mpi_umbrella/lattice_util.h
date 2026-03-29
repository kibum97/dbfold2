#ifndef LATTICE_UTIL_H
#define LATTICE_UTIL_H

#include "globals.h" /* for struct atom */

struct atom;

/* MATRIX_SIZE is set to be 20 in init.h (it's a u. char, not a macro) */

struct cell {
    char         X, Y, Z;
    struct cell *neighbors[27];
    short atom_list[MAX_CELL_ATOMS];  // AB: Each cell has this array called atom_list which has
                                      // MAX_CELL_ATOMS = 100 elements
    unsigned char natoms;
};

void          FindLatticeCoordinates(struct atom *);
void          UpdateLattice(short, short *);
void          InitializeMatrix();
unsigned char PerBound(signed char);
void          CopyLatticeCoordinates(struct atom, struct atom *);

#endif