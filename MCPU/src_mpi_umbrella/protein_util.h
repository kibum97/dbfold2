#ifndef PROTEIN_UTIL_H
#define PROTEIN_UTIL_H

#define SP3_ANGLE 109.5
#define CA_CB_DISTANCE 1.54

Float Phi(struct residue, struct residue, struct atom *);
Float Psi(struct residue, struct residue, struct atom *);
Float CalculateTorsion(struct atom *, int, int, int, int, Float);
void  AddCB(struct atom *, struct residue, struct atom *);
Float RadiusOfGyration(struct atom *protein, int num_atoms);
void  CenterProtein(struct atom **, int);
void  dRms(struct atom *, struct atom *, int, Float *, Float *);
void  SSdRms(struct atom *, struct atom *, int, Float *, Float *, int, int);

#endif