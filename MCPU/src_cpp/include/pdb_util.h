#ifndef PDB_UTIL_H
#define PDB_UTIL_H

#include "globals.h"
#include <string_view>

void ParsePDBLine(
    struct Simulation *sim,
    struct System *sys,
    struct MCIntegrator *integrator,
    struct Topology *top,
    char *line, struct atom *protein, int *Natoms);
int TypeAtom(char *s, char *res);
int IsSidechainAtom(char *atomname);
int IsDesignedResidue(int res);
int IsCoreResidue(int res);
int GetAminoNumber(char *name);
void PrintPDB(
    struct Simulation *sim,
    struct Context *ctx,
    struct Topology *top,
    const char *filename);
void PrintPDB_Emin(
    struct Simulation *sim,     
    struct Context *ctx,
    struct Topology *top,
    const char *filename);
void PrintPDB_RMSDmin(
    struct Simulation *sim,
    struct Context *ctx,
    struct Topology *top,
    const char *filename);
int SMoGType(
    struct Simulation *sim,
    char *atomname, char *residue);
int GetSMoGType(std::string_view c);
void PrintReplica(
    struct Simulation *sim,
    struct Context *ctx,
    struct Topology *top,
    char *filename);
int GetReplica(
    struct Simulation *sim,
    struct Context *ctx,
    struct Topology *top,
    char *filename);

#endif