#ifndef HBONDS_H
#define HBONDS_H

float    HydrogenBonds();
float    FoldHydrogenBonds();
void     InitializeHydrogenBonding();
int      MatchAtomname(char *);
long int CheckHBond(int, int);

#endif