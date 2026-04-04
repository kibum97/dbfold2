#include "contacts.h"

#include "globals.h"
#include "tripep_closure.h"
#include "vector.h"

double contact_calpha_cutoff = 7.; /* CONTACTS DEF*/

void Contacts(struct Context *ctx, const struct Topology *top, const struct System *sys) {
    int i, j;
    /* Re-initializes the contact matrices */
    /* Counts clashes, and native/non-native contacts */
    /* Values for ncontacts, nnon_native_contacts will not be correct */
    /*        until GetNativeContacts has been called once */

    ctx->nclashes  = 0;
    ctx->ncontacts = 0;

    for (i = 0; i < top->natoms; i++) {
        for (j = i + 1; j < top->natoms; j++) {
            ctx->data[i][j].clashes = ctx->data[j][i].clashes = 0;
            ctx->data[i][j].contacts = ctx->data[j][i].contacts = 0;
            if (ctx->data[i][j].check_contacts || ctx->data[i][j].check_clashes)
                CheckForContacts(ctx, sys, i, j);
        }
    }
    return;
}

void TypeContacts(struct Context *ctx, const struct System *sys, const struct Topology *top) {
    int i, j;

    for (i = 0; i < sys->MAX_TYPES; i++)
        for (j = 0; j < sys->MAX_TYPES; j++)
            ctx->type_contacts[i][j] = 0;

    for (i = 0; i < top->natoms; i++)
        for (j = i + 1; j < top->natoms; j++)
            if (ctx->data[i][j].contacts) {
                if (ctx->native[i].smogtype <= ctx->native[j].smogtype)
                    ctx->type_contacts[ctx->native[i].smogtype][ctx->native[j].smogtype]++;
                else
                    ctx->type_contacts[ctx->native[j].smogtype][ctx->native[i].smogtype]++;
            }

    return;
}

void CheckForContacts(struct Context *ctx, const struct System *sys, short a, short b) {
    Vec3i XX = ctx->native[a].xyz_int;
    Vec3i YY = ctx->native[b].xyz_int;
    double distance = (XX - YY).norm();

    if (ctx->data[a][b].check_clashes &&
        distance < sys->hard_core[ctx->native[a].smogtype][ctx->native[b].smogtype]) {
        ctx->data[a][b].clashes = ctx->data[b][a].clashes = 1;
        ctx->nclashes++;
    }

    if (ctx->data[a][b].check_contacts) {
        if ((distance <=
             sys->contact_distance[ctx->native[a].smogtype][ctx->native[b].smogtype].b) &&
            (distance >=
             sys->contact_distance[ctx->native[a].smogtype][ctx->native[b].smogtype].a)) {
            ctx->data[a][b].contacts = ctx->data[b][a].contacts = 1;
            ctx->ncontacts++;
        }
    }
    return;
}

unsigned char CheckForDeltaContacts(struct contact_data *Data, Vec3i XX,
                                    Vec3i YY, short type_a, short type_b,
                                    const struct System *sys, struct monte_carlo_flags mc_flags) {
    double   distance;
    distance = (XX - YY).norm();

    if (Data->check_contacts && Data->check_clashes) {
        if (((distance <= sys->contact_distance[type_a][type_b].b) &&
             (distance >= sys->contact_distance[type_a][type_b].a)) != Data->contacts) {
            Data->delta_contacts = !(Data->contacts);
            if ((distance < sys->hard_core[type_a][type_b]) != (Data->clashes))
                Data->delta_clashes = mc_flags.clashed = !(Data->clashes);
            return 3;
        } else if ((distance < sys->hard_core[type_a][type_b]) != (Data->clashes)) {
            Data->delta_clashes = mc_flags.clashed = !(Data->clashes);
            return 2;
        } else
            return 0;
    } else if (Data->check_clashes) {
        if ((distance < sys->hard_core[type_a][type_b]) != Data->clashes) {
            Data->delta_clashes = mc_flags.clashed = !(Data->clashes);
            return 2;
        } else
            return 0;
    } else if (Data->check_contacts) {
        if (((distance <= sys->contact_distance[type_a][type_b].b) &&
             (distance >= sys->contact_distance[type_a][type_b].a)) != (Data->contacts)) {
            Data->delta_contacts = !(Data->contacts);
            return 1;
        } else
            return 0;
    }

    return 0;
}

//========================================================================================================
//========================================================================================================
void check_bb_contacts(struct Context *ctx, short a, short b, const struct Topology *top,
                       const struct Simulation *sim) {
    Vec3 coo, oxy, cao; // acceptor
    Vec3 nit, con, can, hyd;  // donor
    double        len_ho, ang_h, ang_o, dih_coo;
    Vec3 V3, H;
    int           coo_at, oxy_at, cao_at;
    int           nit_at, con_at, can_at;
    int           passed = 1;

    if ((ctx->native[a].smogtype == top->bb_O_type) &&
        (ctx->native[b].smogtype == top->bb_N_type))  // a == O, b == N
    {
        if (b == 0)
            return;
        coo_at = ctx->native_residue[ctx->native[a].res_num].C;
        oxy_at = ctx->native_residue[ctx->native[a].res_num].O;
        cao_at = ctx->native_residue[ctx->native[a].res_num].CA;
        nit_at = ctx->native_residue[ctx->native[b].res_num].N;
        con_at = ctx->native_residue[ctx->native[b - 1].res_num].C;
        can_at = ctx->native_residue[ctx->native[b].res_num].CA;
    }

    else if ((ctx->native[a].smogtype == top->bb_N_type) &&
             (ctx->native[b].smogtype == top->bb_O_type))  // a ==O, b == N
    {
        if (a == 0)
            return;
        coo_at = ctx->native_residue[ctx->native[b].res_num].C;
        oxy_at = ctx->native_residue[ctx->native[b].res_num].O;
        cao_at = ctx->native_residue[ctx->native[b].res_num].CA;
        nit_at = ctx->native_residue[ctx->native[a].res_num].N;
        con_at = ctx->native_residue[ctx->native[a - 1].res_num].C;
        can_at = ctx->native_residue[ctx->native[a].res_num].CA;
    } else {
        fprintf(sim->STATUS, "Error in bb_interactions!\n");
        exit(1);
    }

    coo = ctx->native[coo_at].xyz;
    oxy = ctx->native[oxy_at].xyz;
    nit = ctx->native[nit_at].xyz;
    cao = ctx->native[cao_at].xyz;
    con = ctx->native[con_at].xyz;
    can = ctx->native[can_at].xyz;

    MakeVector(ctx->native[nit_at].xyz, ctx->native[can_at].xyz, &H);
    MakeVector(ctx->native[nit_at].xyz, ctx->native[con_at].xyz, &V3);
    Add(V3, &H);
    Normalize(&H);
    Inverse(&H);
    Add(ctx->native[nit_at].xyz, &H);
    hyd = H;

    c_bnd_len(hyd, oxy, &len_ho);
    c_bnd_ang(oxy, hyd, nit, &ang_h);
    c_bnd_ang(coo, oxy, hyd, &ang_o);
    c_dih_ang(cao, coo, oxy, hyd, &dih_coo);
    ang_h   = ang_h * rad2deg;
    ang_o   = ang_o * rad2deg;
    dih_coo = dih_coo * rad2deg;

    fprintf(sim->STATUS, "%4d %4d %4d %4d %2d %2d %2d %2d %2d %7.3f %7.3f %7.3f %7.3f\n", a, b,
            oxy_at, nit_at, ctx->native[a].smogtype, ctx->native[b].smogtype,
            ctx->native[a].res_num, ctx->native[b].res_num,
            abs(ctx->native[a].res_num - ctx->native[b].res_num), len_ho, ang_h, ang_o, dih_coo);

    if (passed == 1) {
        ctx->data[a][b].contacts = ctx->data[b][a].contacts = 1;
        ctx->ncontacts++;
    }

    return;
}

void NewDeltaContacts(struct Context *ctx, const struct System *sys, short rotate_natoms,
                      short *rotate_atom, char *not_rotated) {
    int    i, j, k;
    short  temp;
    int    N, O, M, total_pairs = 0;
    short *A = 0;

    for (i = 0; i < rotate_natoms; i++) {
        N               = rotate_atom[i];
        temp_atom       = &ctx->native[N];
        temp            = temp_atom->smogtype;
        temp_xyz_int    = temp_atom->xyz_int;
        temp_cell3      = ctx->prev_native[N].matrix;
        temp_cell_array = temp_atom->matrix->neighbors;
        for (j = 0; j < 27; j++) { /* loops over new neighbors and determines changes */
            temp_cell = temp_cell_array[j];
            if (temp_cell->natoms) {
                k = 0;
                O = temp_cell->natoms;
                A = temp_cell->atom_list;
                while (k < O) { /* this is O not zero */
                    M = A[k];
                    if (is_rotated[N] != is_rotated[M] &&
                        (ctx->data[N][M].check_contacts || ctx->data[N][M].check_clashes)) {
                        if (CheckForDeltaContacts(&ctx->data[N][M], temp_xyz_int,
                                                  ctx->native[M].xyz_int, temp,
                                                  ctx->native[M].smogtype, sys, ctx->mc_flags)) {
                            ctx->ab[total_pairs].a   = N;
                            ctx->ab[total_pairs++].b = M;
                        }
                        if (ctx->mc_flags.clashed)
                            return;
                    }
                    k++;
                }
            }
        }  // for (j=0; j<27; j++)

        /* if atom has moved out of its cell, loop over old neighbors as well... */
        if (temp_cell_array[13] != temp_cell3->neighbors[13])
            for (j = 0; j < 27; j++) {
                k = 0;
                while (k < 27 && (temp_cell_array[k] != temp_cell3->neighbors[j]))
                    k++;
                if (k == 27) { /* if an old neighbor cell is no longer a neighbor cell... */
                    temp_cell = temp_cell3->neighbors[j];
                    if (temp_cell->natoms) {
                        k = 0;
                        O = temp_cell->natoms;
                        A = temp_cell->atom_list;
                        while (k < O) {
                            M = A[k];
                            /* if an atom in such a cell used to be in contact with N... */
                            if (is_rotated[N] != is_rotated[M] && ctx->data[N][M].check_contacts &&
                                ctx->data[N][M].contacts) {
                                /* ... then it is no longer in contact! */
                                /* so add pair to list, and default delta_contact = 0 turns off the
                                 * contact */
                                ctx->ab[total_pairs].a   = N;
                                ctx->ab[total_pairs++].b = M;
                            }
                            k++;
                        }
                    }
                }
            }  // for (j=0; j<27; j++)

    }  // for (i=0; i<rotated_natoms; i++)

    return;
}

////////////////////////////////////////////////////////
/* wmj */
///////////////////////////////////////////////////////////

// const double contact_calpha_cutoff = 7.; /* CONTACTS DEF*/
// AB commented out the above and made it a global in backbone.h, so that it's 7 by default but can
// be changed
void fill_calpha_contact_string(const struct Topology *top, const struct residue *residues,
                                const struct atom *atoms, short *contactstring) {
    int    i, j, atomi, atomj;
    double xi, yi, zi, xj, yj, zj, dist;
    for (i = 0; i < top->nresidues; i++) {
        atomi = residues[i].CA;
        for (j = i + top->min_seq_sep; j < top->nresidues; j++) /* CONTACTS DEF*/
        {
            atomj = residues[j].CA;
            dist  = (atoms[atomi].xyz - atoms[atomj].xyz).norm();
            if (dist < contact_calpha_cutoff * contact_calpha_cutoff) {
                contactstring[i * top->nresidues + j] = contactstring[j * top->nresidues + i] = 1;
            } else {
                contactstring[i * top->nresidues + j] = contactstring[j * top->nresidues + i] = 0;
            }
        }
    }
}

double number_of_calpha_contacts(const struct Topology *top, const struct residue *residues,
                                 const struct atom *atoms) {
    int    i, j, atomi, atomj;
    double xi, yi, zi, xj, yj, zj, dist;
    double number_of_contacts = 0.;
    for (i = 0; i < top->nresidues; i++) {
        atomi = residues[i].CA;
        for (j = i + top->min_seq_sep; j < top->nresidues; j++) /* CONTACTS DEF*/
        {
            atomj = residues[j].CA;
            dist  = (atoms[atomi].xyz - atoms[atomj].xyz).norm();
            if (dist < contact_calpha_cutoff * contact_calpha_cutoff) {
                number_of_contacts += 1.;
            }
        }
    }
    return number_of_contacts;
}

double number_of_calpha_native_contacts(const struct Topology *top, const struct residue *residues,
                                        const struct atom *atoms, const short *contactstring) {
    int    i, j, atomi, atomj;
    double xi, yi, zi, xj, yj, zj, dist;
    double number_of_contacts = 0.;
    for (i = 0; i < top->nresidues; i++) {
        atomi = residues[i].CA;
        for (j = i + top->min_seq_sep; j < top->nresidues; j++) /* CONTACTS DEF*/
        {
            atomj = residues[j].CA;
            dist  = (atoms[atomi].xyz - atoms[atomj].xyz).norm();
            if (dist < contact_calpha_cutoff * contact_calpha_cutoff &&
                contactstring[top->nresidues * i + j]) {
                number_of_contacts += 1.;
            }
        }
    }
    return number_of_contacts;
}

double number_of_calpha_nonnative_contacts(const struct Topology *top,
                                           const struct residue *residues, const struct atom *atoms,
                                           const short *contactstring) {
    int    i, j, atomi, atomj;
    double xi, yi, zi, xj, yj, zj, dist;
    double number_of_contacts = 0.;
    for (i = 0; i < top->nresidues; i++) {
        atomi = residues[i].CA;
        for (j = i + top->min_seq_sep; j < top->nresidues; j++) /* CONTACTS DEF*/
        {
            atomj = residues[j].CA;
            dist  = (atoms[atomi].xyz - atoms[atomj].xyz).norm();
            if (dist < contact_calpha_cutoff * contact_calpha_cutoff &&
                !contactstring[top->nresidues * i + j]) {
                number_of_contacts += 1.;
            }
        }
    }
    return number_of_contacts;
}

double hamming_distance_calpha_contacts(const struct Topology *top, const struct residue *residues,
                                        const struct atom *atoms, const short *contactstring) {
    int    i, j, atomi, atomj;
    double xi, yi, zi, xj, yj, zj, dist;
    short  iscontact;
    double hamming_distance = 0.;
    for (i = 0; i < top->nresidues; i++) {
        atomi = residues[i].CA;
        for (j = i + top->min_seq_sep; j < top->nresidues; j++) /* CONTACTS DEF*/
        {
            atomj     = residues[j].CA;
            dist      = (atoms[atomi].xyz - atoms[atomj].xyz).norm();
            iscontact = dist < contact_calpha_cutoff * contact_calpha_cutoff;
            hamming_distance += iscontact ^ contactstring[top->nresidues * i + j];
        }
    }
    return hamming_distance;
}
/* end wmj */