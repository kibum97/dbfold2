#ifndef PDB_UTIL_H
#define PDB_UTIL_H

#include "globals.h"

void ParsePDBLine(char *, struct atom *, int *);
int  TypeAtom(char *, char *);
int  IsSidechainAtom(char *);
int  IsDesignedResidue(int);
int  IsCoreResidue(int);
int  GetAminoNumber(char *);
void PrintPDB(char *);
void PrintPDB_Emin(char *);
void PrintPDB_RMSDmin(char *);
int  SMoGType(char *, char *);
int  GetSMoGType(char *);
void PrintReplica(char *);
int  GetReplica(char *);

#endif