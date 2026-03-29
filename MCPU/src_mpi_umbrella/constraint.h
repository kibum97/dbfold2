#ifndef CONSTRAINT_H
#define CONSTRAINT_H

struct alignment constraint_align;

void  Init_frag_rmsd_constraints();
float Compute_constraint_energy(struct residue *res, struct atom *at);
float Compute_frag_rmsd_constraint_energy(struct backbone  struct1[MAXSEQUENCE],
                                          struct backbone  struct2[MAXSEQUENCE],
                                          struct alignment algn);
int   Argmin(float distances[]);
void  Read_constraints();

#endif