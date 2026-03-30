#ifndef INIT_H
#define INIT_H

#include "atom.h"
#include "globals.h"

void SetProgramOptions(int argc, char *argv[]);
void InitializeProtein();
void DetermineTriplets();
void SetRadii();
void SetHardCore();
void SetContactDistance();
void ReadNative(char *file_name, struct atom *protein, int *Natoms);
void GetResidueInfo(struct atom *, struct residue *, int, int);
void GetPhiPsi(struct atom *, struct residue *, int);
void GetChi();
void InitializeData();
void ReadSidechainTorsionData();
void InitializeSidechainRotationData();
void InitializeBackboneRotationData();
void CheckCorrelation(struct contact_data **, struct atom *, struct residue *, int);
void ReadPotential();
void ReadDistPotential();
int  SkipSelf(int, int, struct atom *, struct residue *);
int  SkipNeighbors(int, int, struct atom *, struct residue *);
void TurnOffNativeClashes(int);
void ReadAlignment();
void SetupAlignmentStructure();
void SetupAlignmentPotential();
void ReadAvgChis(void);
void ReadHelicityData();
void InitializeHydrogenBonding();
int  MatchAtomname(char *);

#endif
