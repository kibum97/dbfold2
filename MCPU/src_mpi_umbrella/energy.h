#ifndef ENERGY_H
#define ENERGY_H

Float FullAtomEnergy();
void  ResetEnergies(long int check);
Float HydrogenBonds();
Float torsionenergy();
Float sctenergy();
void  aromatic_center(int, struct vector *);
float aromatic_plane(int, int, struct vector *, struct vector *);
Float aromaticenergy();

Float Rgenergy();

#endif