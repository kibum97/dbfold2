#ifndef MOVE_H
#define MOVE_H

#include "globals.h" /* for global variables */

void BackboneMove(struct Context *ctx, struct Topology *top, struct MCIntegrator *integrator,
                  float step_size);
void LoopBackboneMove(struct Context *ctx, struct Simulation *sim,
                      struct MCIntegrator *integrator, struct Topology *top,
                      struct System *sys, float absolute_step_size);
void LocalBackboneMove(struct Context *ctx, struct Topology *top, struct System *sys,
                       struct Simulation *sim, struct MCIntegrator *integrator, float step_size);
void SidechainMove(struct Context *ctx, struct Topology *top, struct MCIntegrator *integrator,
                   struct System *sys);
void MakeSidechainMove(struct Context *ctx, struct Topology *top,
                       struct MCIntegrator *integrator, struct System *sys);
void MakeMove(struct Context *ctx, struct System *sys, struct Topology *top,
              struct Simulation *sim, struct MCIntegrator *integrator, float step_size,
              float use_global_bb_moves);

#endif