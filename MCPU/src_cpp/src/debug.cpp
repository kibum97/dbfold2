#include "debug.h"

#include <stdio.h>

#if DEBUG
void Debug(
    struct Context *ctx,
    struct Topology *top,
    struct Simulation *sim,
    struct System *sys
) {
    int   i, j, k;
    short temp;

    DebugContacts(sim, top, ctx, sys);
    for (i = 0; i < top->natoms; i++)
        for (j = i + 1; j < top->natoms; j++)
            if (ctx->data[i][j].clashes != sim->debug_clashes[i][j])
                fprintf(sim->STATUS, "CLASH_MISMATCH\t(%d, %d)\t%d\t%d\n", i, j, ctx->data[i][j].clashes,
                        sim->debug_clashes[i][j]);
    for (i = 0; i < top->natoms; i++)
        for (j = i + 1; j < top->natoms; j++)
            if (ctx->data[i][j].contacts != sim->debug_contacts[i][j]) {
                fprintf(sim->STATUS, "CONTACT_MISMATCH\t(%d, %d)\t%d\t%d", i, j, ctx->data[i][j].contacts,
                        sim->debug_contacts[i][j]);
                for (k = 0; k < ctx->total_pairs; k++)
                    if (ctx->ab[k].a == i && ctx->ab[k].b == j) {
                        fprintf(sim->STATUS, "%d\n", k);
                        break;
                    }
                if (k == ctx->total_pairs)
                    fprintf(sim->STATUS, "ERROR\n");
            }
    temp = 0;
    for (i = 0; i < top->natoms; i++)
        for (j = i + 1; j < top->natoms; j++)
            temp += ctx->data[i][j].contacts;
    fprintf(sim->STATUS, "REPORTED\t%d\t%d\t%d\n", ctx->nclashes, ctx->ncontacts, temp);
    temp = 0;
    for (i = 0; i < top->natoms; i++)
        for (j = i + 1; j < top->natoms; j++)
            temp += sim->debug_contacts[i][j];
    fprintf(sim->STATUS, "DEBUG\t\t%d\t%d\t%d\n", sim->debug_nclashes, sim->debug_ncontacts, temp);

    return;
}

void DebugContacts(
    struct Simulation *sim,
    struct Topology *top,
    struct Context *ctx,
    struct System *sys
) {
    int i, j;
    sim->debug_nclashes  = 0;
    sim->debug_ncontacts = 0;
    for (i = 0; i < top->natoms; i++)
        for (j = i + 1; j < top->natoms; j++) {
            sim->debug_clashes[i][j]  = 0;
            sim->debug_contacts[i][j] = 0;
            if (ctx->data[i][j].check_contacts || ctx->data[i][j].check_clashes)
                CheckForDebugContacts(sim, top, ctx, sys, i, j);
        }

    return;
}

void CheckForDebugContacts(
    struct Simulation *sim,
    struct Topology *top,
    struct Context *ctx,
    struct System *sys,
    int a, int b) {
    distance =
        (ctx->native[a].xyz_int.x - ctx->native[b].xyz_int.x) * (ctx->native[a].xyz_int.x - ctx->native[b].xyz_int.x) +
        (ctx->native[a].xyz_int.y - ctx->native[b].xyz_int.y) * (ctx->native[a].xyz_int.y - ctx->native[b].xyz_int.y) +
        (ctx->native[a].xyz_int.z - ctx->native[b].xyz_int.z) * (ctx->native[a].xyz_int.z - ctx->native[b].xyz_int.z);

    if (ctx->data[a][b].check_clashes && distance < sys->hard_core[ctx->native[a].smogtype][ctx->native[b].smogtype]) {
        sim->debug_clashes[a][b] = sim->debug_clashes[b][a] = 1;
        sim->debug_nclashes++;
    }
    if (ctx->data[a][b].check_contacts &&
        (distance <= sys->contact_distance[ctx->native[a].smogtype][ctx->native[b].smogtype].b) &&
        (distance >= sys->contact_distance[ctx->native[a].smogtype][ctx->native[b].smogtype].a)) {
        sim->debug_contacts[a][b] = sim->debug_contacts[b][a] = 1;
        sim->debug_ncontacts++;
    }

    if (sim->debug_contacts[a][b] != ctx->data[a][b].contacts)
        fprintf(sim->STATUS, "CONTACTS %d, %d --> %f\t%f\t%f\n", a, b, sqrt(distance) / 100,
                sqrt(sys->hard_core[ctx->native[a].smogtype][ctx->native[b].smogtype]) / 100,
                sqrt(sys->contact_distance[ctx->native[a].smogtype][ctx->native[b].smogtype].b) / 100);
    if (sim->debug_clashes[a][b] != ctx->data[a][b].clashes)
        fprintf(sim->STATUS, "CLASHES %d, %d --> %ld\t%ld, %ld\n", a, b, distance,
                sys->hard_core[ctx->native[a].smogtype][ctx->native[b].smogtype],
                sys->contact_distance[ctx->native[a].smogtype][ctx->native[b].smogtype].b);

    return;
}

#endif
