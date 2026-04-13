#ifndef INIT_H
#define INIT_H

#include "atom.h"
#include "globals.h"

void SetProgramOptions(struct Simulation *sim, struct System *sys, struct MCIntegrator *integrator,
                       struct Context *ctx, struct Topology *top, int argc, char *argv[]);
void InitializeProtein(struct Context *ctx, struct Topology *top, struct Simulation *sim,
                       struct System *sys, struct MCIntegrator *integrator);
void DetermineTriplets(struct Topology *top, struct MCIntegrator *integrator,
                       struct Simulation *sim, struct System *sys);
void SetRadii(struct System *sys);
void SetHardCore(struct System *sys, struct Topology *top, struct Context *ctx);
void SetContactDistance(struct System *sys, struct Topology *top, struct Context *ctx);
void SetupMuPotential(struct System *sys, struct Topology *top, struct Context *ctx,
                      struct Simulation *sim);
void ReadNative(struct Simulation *sim,
    struct System *sys, struct MCIntegrator *integrator, struct Topology *top,
    const char *file_name, struct atom *protein, int *Natoms);
void ReadTypesFile(struct Simulation *sim, struct Topology *top, struct System *sys);
void GetResidueInfo(struct atom *Chain, struct residue *Residue, int Nres, int Natoms);
void GetPhiPsi(struct atom *Chain, struct residue *Residue, int Nres);
void GetChi(struct Topology *top, struct Context *ctx, struct MCIntegrator *integrator);
void InitializeData(struct Context *ctx, struct Topology *top, struct Simulation *sim,
                    struct System *sys);
void ReadSidechainTorsionData(struct Simulation *sim, struct Topology *top);
void InitializeSidechainRotationData(struct Topology *top, struct MCIntegrator *integrator,
                                     struct Context *ctx, struct System *sys,
                                     struct Simulation *sim);
void InitializeBackboneRotationData(struct Topology *top, struct MCIntegrator *integrator,
                                    struct Context *ctx);
void CheckCorrelation(struct contact_data **Data, struct atom *Protein, struct residue *Residue,
                      int Natoms, struct Simulation *sim);
void ReadPotential(struct Simulation *sim, struct System *sys);
// void ReadDistPotential();
int  SkipSelf(int s, int b, struct atom *Protein, struct residue *Residue);
int  SkipNeighbors(int i, int j, struct atom *Protein, struct residue *Residue);
void TurnOffNativeClashes(struct Context *ctx, struct Topology *top, struct Simulation *sim,
                          struct System *sys, int ReportBack);
void ReadAvgChis(struct Context *ctx, struct Topology *top, struct Simulation *sim,
                 struct MCIntegrator *integrator, struct System *sys);

#endif
