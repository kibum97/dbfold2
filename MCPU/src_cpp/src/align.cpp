#include "align.h"

#include <string.h>
#include <fstream>

#include "init.h"
#include "pdb_util.h"
#include "vector.h"

void SetAlignHardCore(
    struct System *sys,
    struct Topology *top,
    struct Context *ctx
) {
    int   i, j, k;
    float temp;
    float rad1, rad2;

    if (!sys->USE_GO_POTENTIAL) {
        sys->align_hard_core = (float **)calloc(sys->MAX_TYPES, sizeof(float *));
        for (i = 0; i < sys->MAX_TYPES; i++)
            sys->align_hard_core[i] = (float *)calloc(sys->MAX_TYPES, sizeof(float));

        for (i = 0; i < sys->MAX_TYPES; i++)
            for (j = 0; j < sys->MAX_TYPES; j++) {
                for (k = 0; k < top->natom_type_list; ++k)
                    if (i == top->atom_type_list[k].type_num)
                        break;
                rad1 = sys->radii[TypeAtom(top->atom_type_list[k].atom_name, top->atom_type_list[k].res_name)];
                for (k = 0; k < top->natom_type_list; ++k)
                    if (j == top->atom_type_list[k].type_num)
                        break;
                rad2 = sys->radii[TypeAtom(top->atom_type_list[k].atom_name, top->atom_type_list[k].res_name)];

                temp                  = sys->ALPHA * (rad1 + rad2);
                sys->align_hard_core[i][j] = temp * temp;
            }
    } else {
        sys->align_hard_core = (float **)calloc(top->struct_natoms, sizeof(float *));
        for (i = 0; i < top->struct_natoms; i++)
            sys->align_hard_core[i] = (float *)calloc(top->struct_natoms, sizeof(float));

        for (i = 0; i < top->struct_natoms; i++)
            for (j = 0; j < top->struct_natoms; j++) {
                rad1 = sys->radii[ctx->struct_native[i].atomtype];
                rad2 = sys->radii[ctx->struct_native[j].atomtype];

                temp                  = sys->ALPHA * (rad1 + rad2);
                sys->align_hard_core[i][j] = temp * temp;
            }
    }

    return;
}

void SetAlignContactDistance(
    struct System *sys,
    struct Topology *top,
    struct Context *ctx
) {
    float rad1, rad2;
    float temp;
    int   i, j, k;

    if (!sys->USE_GO_POTENTIAL) {
        sys->align_con_dist = (struct align_cutoff **)calloc(sys->MAX_TYPES, sizeof(struct align_cutoff *));
        for (i = 0; i < sys->MAX_TYPES; i++)
            sys->align_con_dist[i] =
                (struct align_cutoff *)calloc(sys->MAX_TYPES, sizeof(struct align_cutoff));

        for (i = 0; i < sys->MAX_TYPES; i++)
            for (j = 0; j < sys->MAX_TYPES; j++) {
                for (k = 0; k < top->natom_type_list; ++k)
                    if (i == top->atom_type_list[k].type_num)
                        break;
                rad1 = sys->radii[TypeAtom(top->atom_type_list[k].atom_name, top->atom_type_list[k].res_name)];
                for (k = 0; k < top->natom_type_list; ++k)
                    if (j == top->atom_type_list[k].type_num)
                        break;
                rad2 = sys->radii[TypeAtom(top->atom_type_list[k].atom_name, top->atom_type_list[k].res_name)];

                temp                   = sys->LAMBDA * sys->ALPHA * (rad1 + rad2);
                sys->align_con_dist[i][j].b = temp * temp;
                sys->align_con_dist[i][j].a = 0;
            }

        if (sys->MAX_TYPES == 84) {
            sys->align_con_dist[79][82].b = sys->align_con_dist[82][79].b = 3.25 * 3.25;
            sys->align_con_dist[79][83].b = sys->align_con_dist[83][79].b = 3.25 * 3.25;
            sys->align_con_dist[79][82].a = sys->align_con_dist[82][79].a = 2.75 * 2.75;
            sys->align_con_dist[79][83].a = sys->align_con_dist[83][79].a = 2.75 * 2.75;
            for (j = 79; j <= 83; ++j)
                for (i = 0; i < sys->MAX_TYPES; ++i)
                    if ((i < 79) || (i > 83))
                        sys->align_con_dist[i][j].b = sys->align_con_dist[j][i].b = 0;
        }
        if (sys->MAX_TYPES == 28) {
            sys->align_con_dist[23][26].b = sys->align_con_dist[26][23].b = 3.25 * 3.25;
            sys->align_con_dist[23][27].b = sys->align_con_dist[27][23].b = 3.25 * 3.25;
            sys->align_con_dist[23][26].a = sys->align_con_dist[26][23].a = 2.75 * 2.75;
            sys->align_con_dist[23][27].a = sys->align_con_dist[27][23].a = 2.75 * 2.75;
            for (j = 23; j <= 27; ++j)
                for (i = 0; i < sys->MAX_TYPES; ++i)
                    if ((i < 23) || (i > 27))
                        sys->align_con_dist[i][j].b = sys->align_con_dist[j][i].b = 0;
        }
    } else {
        sys->align_con_dist =
            (struct align_cutoff **)calloc(top->struct_natoms, sizeof(struct align_cutoff *));
        for (i = 0; i < top->struct_natoms; i++)
            sys->align_con_dist[i] =
                (struct align_cutoff *)calloc(top->struct_natoms, sizeof(struct align_cutoff));

        for (i = 0; i < top->struct_natoms; i++)
            for (j = 0; j < top->struct_natoms; j++) {
                rad1 = sys->radii[ctx->struct_native[i].atomtype];
                rad2 = sys->radii[ctx->struct_native[j].atomtype];

                temp                   = sys->LAMBDA * sys->ALPHA * (rad1 + rad2);
                sys->align_con_dist[i][j].b = temp * temp;
                sys->align_con_dist[i][j].a = 0;
            }
    }

    return;
}

void ReadAlignment(
    struct Topology *top,
    struct System *sys
) {
    int seq1, seq2, struct1, struct2, i;
    int nalign_seq, nalign_struct;

    top->map_to_seq    = (int *)calloc(1000, sizeof(int));
    top->map_to_struct = (int *)calloc(1000, sizeof(int));
    nalign_seq = nalign_struct = 0;
    sys->nseg                       = 0;
    seq1                       = 0;
    seq2                       = top->nresidues - 1;
    struct1                    = 0;
    struct2                    = top->nresidues - 1;
    for (i = seq1; i <= seq2; ++i)
        top->map_to_seq[nalign_seq++] = i;
    for (i = struct1; i <= struct2; ++i)
        top->map_to_struct[nalign_struct++] = i;
    sys->str_segment[sys->nseg].a = struct1;
    sys->str_segment[sys->nseg].b = struct2;
    sys->seq_segment[sys->nseg].a = seq1;
    sys->seq_segment[sys->nseg].b = seq2;
    ++sys->nseg;
    top->nalign = nalign_seq;
}

void AlignCheckForContacts(
    struct Context *ctx,
    struct System *sys,
    short a, short b) {
    float distance;

    distance = D2(ctx->struct_native[a].xyz, ctx->struct_native[b].xyz);

    if (ctx->struct_data[a][b].check_clashes &&
        distance < sys->align_hard_core[ctx->struct_native[a].smogtype][ctx->struct_native[b].smogtype]) {
        ctx->struct_data[a][b].clashes = ctx->struct_data[b][a].clashes = 1;
        ctx->struct_nclashes++;
    }

    if (ctx->struct_data[a][b].check_contacts)
        if ((distance <= sys->align_con_dist[ctx->struct_native[a].smogtype][ctx->struct_native[b].smogtype].b) &&
            (distance >= sys->align_con_dist[ctx->struct_native[a].smogtype][ctx->struct_native[b].smogtype].a)) {
            ctx->struct_data[a][b].contacts = ctx->struct_data[b][a].contacts = 1;
            ctx->struct_ncontacts++;
        }

    return;
}

void AlignContacts(
    struct Context *ctx,
    struct Topology *top,
    struct System *sys
) {
    int i, j;
    /* Re-initializes the contact matrices */
    /* Counts clashes, and native/non-native contacts */
    /* Values for ncontacts, nnon_native_contacts will not be correct */
    /*        until GetNativeContacts has been called once */

    ctx->struct_nclashes  = 0;
    ctx->struct_ncontacts = 0;

    for (i = 0; i < top->struct_natoms; i++) {
        for (j = i + 1; j < top->struct_natoms; j++) {
            ctx->struct_data[i][j].clashes = ctx->struct_data[j][i].clashes = 0;
            ctx->struct_data[i][j].contacts = ctx->struct_data[j][i].contacts = 0;
            if (ctx->struct_data[i][j].check_contacts || ctx->struct_data[i][j].check_clashes)
                AlignCheckForContacts(ctx, sys, i, j);
        }
    }
    return;
}

void SetupAlignmentStructure(
    struct Context *ctx,
    struct Topology *top,
    struct System *sys,
    struct MCIntegrator *integrator,
    struct Simulation *sim
) {
    int i;

    ctx->struct_native = (struct atom *)calloc(MAX_ATOMS, sizeof(struct atom));
    ReadNative(sim, sys, integrator, top, sim->structure_file.c_str(), ctx->struct_native, &top->struct_natoms);
    /* despite name, struct_natoms is an int */

    /* this is so unclear:
     * both functions below use struct_native struct_natoms and a bunch of globals
     */
    SetAlignHardCore(sys, top, ctx);
    SetAlignContactDistance(sys, top, ctx);

    top->struct_nresidues = 0;
    for (i = 0; i < top->struct_natoms; i++)
        if (!strncmp(ctx->struct_native[i].atomname, "CA", 2))
            top->struct_nresidues++;
    ctx->struct_residue = (struct residue *)calloc(top->struct_nresidues, sizeof(struct residue));
    GetResidueInfo(ctx->struct_native, ctx->struct_residue, top->struct_nresidues, top->struct_natoms);
    ctx->struct_data = (struct contact_data **)calloc(top->struct_natoms, sizeof(struct contact_data *));
    for (i = 0; i < top->struct_natoms; i++)
        ctx->struct_data[i] = (struct contact_data *)calloc(top->struct_natoms, sizeof(struct contact_data));

    GetPhiPsi(ctx->struct_native, ctx->struct_residue, top->struct_nresidues);
    CheckCorrelation(ctx->struct_data, ctx->struct_native, ctx->struct_residue, top->struct_natoms, sim);

    top->seq_to_struct = (int *)calloc(top->struct_nresidues, sizeof(int));
    top->struct_to_seq = (int *)calloc(top->struct_nresidues, sizeof(int));

    for (i = 0; i < top->struct_nresidues; ++i)
        top->seq_to_struct[i] = -1;
    for (i = 0; i < top->struct_nresidues; ++i)
        top->struct_to_seq[i] = -1;

    for (i = 0; i < top->nalign; ++i) {
        top->seq_to_struct[top->map_to_seq[i]]    = top->map_to_struct[i];
        top->struct_to_seq[top->map_to_struct[i]] = top->map_to_seq[i];
    }

    AlignContacts(ctx, top, sys);
    fprintf(sim->STATUS, "---TEMPLATE---\n");
    fprintf(sim->STATUS, "  # of clashes:\t\t%d\n  # of contacts:\t%d\n\n", ctx->struct_nclashes,
            ctx->struct_ncontacts);
}

void SetupAlignmentPotential(
    struct Topology *top,
    struct Context *ctx,
    struct System *sys,
    struct Simulation *sim
) {
    int            ngo_seq, ngo_str;
    int            i, j, k, l, m, str_bbA[4], str_bbB[4], seq_bbA[4], seq_bbB[4];
    int          **done_pairs;
    struct residue str_resA, str_resB, seq_resA, seq_resB;
    int            strA, strB, seqA, seqB;
    short          same_segment;

    /* setup aligned backbone potential */
    ngo_seq = 0;
    ngo_str = 0;
    for (i = 0; i < top->nalign; ++i)
        for (j = i + 1; j < top->nalign; ++j) {
            str_resA   = ctx->struct_residue[top->map_to_struct[i]];
            str_resB   = ctx->struct_residue[top->map_to_struct[j]];
            seq_resA   = ctx->native_residue[top->map_to_seq[i]];
            seq_resB   = ctx->native_residue[top->map_to_seq[j]];
            str_bbA[0] = str_resA.CA;
            str_bbB[0] = str_resB.CA;
            str_bbA[1] = str_resA.N;
            str_bbB[1] = str_resB.N;
            str_bbA[2] = str_resA.C;
            str_bbB[2] = str_resB.C;
            str_bbA[3] = str_resA.O;
            str_bbB[3] = str_resB.O;
            seq_bbA[0] = seq_resA.CA;
            seq_bbB[0] = seq_resB.CA;
            seq_bbA[1] = seq_resA.N;
            seq_bbB[1] = seq_resB.N;
            seq_bbA[2] = seq_resA.C;
            seq_bbB[2] = seq_resB.C;
            seq_bbA[3] = seq_resA.O;
            seq_bbB[3] = seq_resB.O;

            /* determine whether i and j are in same aligned segment in structure */

            same_segment = 0;
            for (m = 0; m < sys->nseg; ++m)
                if ((top->map_to_struct[i] >= sys->str_segment[m].a) &&
                    (top->map_to_struct[i] <= sys->str_segment[m].b) &&
                    (top->map_to_struct[j] >= sys->str_segment[m].a) &&
                    (top->map_to_struct[j] <= sys->str_segment[m].b))
                    same_segment = 1;

            /* generate potential: */
            /*   if aligned residues are in contact, assign attraction  */
            /*   if not, then assign repulsion only if they are in same segment */
            /*   otherwise, assign non specific energy */

            for (k = 0; k < 4; ++k)
                for (l = 0; l < 4; ++l)
                    if (ctx->struct_data[str_bbA[k]][str_bbB[l]].contacts) {
                        sys->potential[seq_bbA[k]][seq_bbB[l]] = sys->NATIVE_ATTRACTION;
                        sys->potential[seq_bbB[l]][seq_bbA[k]] = sys->NATIVE_ATTRACTION;
                        ++ngo_str;
                        ++ngo_seq;
                    } else if (same_segment) {
                        sys->potential[seq_bbA[k]][seq_bbB[l]] = sys->NON_NATIVE_REPULSION;
                        sys->potential[seq_bbB[l]][seq_bbA[k]] = sys->NON_NATIVE_REPULSION;
                    } else {
                        sys->potential[seq_bbA[k]][seq_bbB[l]] = sys->NON_SPECIFIC_ENERGY;
                        sys->potential[seq_bbB[l]][seq_bbA[k]] = sys->NON_SPECIFIC_ENERGY;
                    }
        }
    fprintf(sim->STATUS, "Sequence Go, backbone: %d\n", ngo_seq);
    fprintf(sim->STATUS, "Structure Go, backbone: %d\n", ngo_str);

    /* setup alignment side-chain potential */

    done_pairs = (int **)calloc(top->nresidues, sizeof(int *));
    for (i = 0; i < top->nresidues; ++i)
        done_pairs[i] = (int *)calloc(top->nresidues, sizeof(int));

    for (i = 0; i < top->struct_natoms; ++i)
        for (j = i + 1; j < top->struct_natoms; ++j) {
            strA = ctx->struct_native[i].res_num;
            strB = ctx->struct_native[j].res_num;
            seqA = top->struct_to_seq[strA];
            seqB = top->struct_to_seq[strB];
            if ((seqA != -1) && (seqB != -1))
                if (ctx->struct_native[i].is_sidechain && ctx->struct_native[j].is_sidechain) {
                    same_segment = 0;
                    for (m = 0; m < sys->nseg; ++m)
                        if ((strA >= sys->str_segment[m].a) && (strA <= sys->str_segment[m].b) &&
                            (strB >= sys->str_segment[m].a) && (strB <= sys->str_segment[m].b))
                            same_segment = 1;
                    if (ctx->struct_data[i][j].contacts) {
                        ++ngo_str;
                        if (!done_pairs[seqA][seqB]) {
                            done_pairs[seqA][seqB] = 1;
                            done_pairs[seqB][seqA] = 1;
                            for (k = 0; k < top->nresidues; ++k)
                                for (l = k + 1; l < top->nresidues; ++l)
                                    if ((ctx->native[k].res_num == seqA) && (ctx->native[l].res_num == seqB))
                                        if (ctx->native[k].is_sidechain && ctx->native[l].is_sidechain) {
                                            sys->potential[k][l] = sys->NATIVE_ATTRACTION;
                                            sys->potential[l][k] = sys->NATIVE_ATTRACTION;
                                            ++ngo_seq;
                                        }
                        }
                    } else {
                        if ((!done_pairs[seqA][seqB]) && (same_segment)) {
                            done_pairs[seqA][seqB] = 1;
                            done_pairs[seqB][seqA] = 1;
                            for (k = 0; k < top->nresidues; ++k)
                                for (l = k + 1; l < top->nresidues; ++l)
                                    if ((ctx->native[k].res_num == seqA) && (ctx->native[l].res_num == seqB))
                                        if (ctx->native[k].is_sidechain && ctx->native[l].is_sidechain) {
                                            sys->potential[k][l] = sys->NON_NATIVE_REPULSION;
                                            sys->potential[l][k] = sys->NON_NATIVE_REPULSION;
                                        }
                        }
                    }
                }
        }

    fprintf(sim->STATUS, "Sequence Go, all: %d\n", ngo_seq);
    fprintf(sim->STATUS, "Structure Go, all: %d\n", ngo_str);

    for (i = 0; i < top->nresidues; ++i)
        free(done_pairs[i]);
    free(done_pairs);
}
