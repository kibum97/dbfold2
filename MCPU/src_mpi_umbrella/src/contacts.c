#include "contacts.h"

#include "tripep_closure.h"
#include "vector.h"

double contact_calpha_cutoff = 7.; /* CONTACTS DEF*/

void Contacts() {
    int i, j;
    /* Re-initializes the contact matrices */
    /* Counts clashes, and native/non-native contacts */
    /* Values for ncontacts, nnon_native_contacts will not be correct */
    /*        until GetNativeContacts has been called once */

    nclashes  = 0;
    ncontacts = 0;

    for (i = 0; i < natoms; i++) {
        for (j = i + 1; j < natoms; j++) {
            data[i][j].clashes = data[j][i].clashes = 0;
            data[i][j].contacts = data[j][i].contacts = 0;
            if (data[i][j].check_contacts || data[i][j].check_clashes)
                CheckForContacts(i, j);
        }
    }
    return;
}

void TypeContacts() {
    int i, j;

    for (i = 0; i < MAX_TYPES; i++)
        for (j = 0; j < MAX_TYPES; j++)
            type_contacts[i][j] = 0;

    for (i = 0; i < natoms; i++)
        for (j = i + 1; j < natoms; j++)
            if (data[i][j].contacts) {
                if (native[i].smogtype <= native[j].smogtype)
                    type_contacts[native[i].smogtype][native[j].smogtype]++;
                else
                    type_contacts[native[j].smogtype][native[i].smogtype]++;
                //	if(((native[i].smogtype==bb_O_type)&&(native[j].smogtype==bb_N_type))||
                //	   ((native[i].smogtype==bb_N_type)&&(native[j].smogtype==bb_O_type)))
                //          fprintf(STATUS, "%4d %4d %2d %2d %2d %2d %2d %3d\n",
                //	      i, j, native[i].smogtype, native[j].smogtype, native[i].res_num,
                // native[j].res_num, 	      abs(native[i].res_num-native[j].res_num),
                // type_contacts[native[j].smogtype][native[i].smogtype]);
            }

    return;
}

void CheckForContacts(short a, short b) {
    double distance =
        (native[a].xyz_int.x - native[b].xyz_int.x) * (native[a].xyz_int.x - native[b].xyz_int.x) +
        (native[a].xyz_int.y - native[b].xyz_int.y) * (native[a].xyz_int.y - native[b].xyz_int.y) +
        (native[a].xyz_int.z - native[b].xyz_int.z) * (native[a].xyz_int.z - native[b].xyz_int.z);
    //  fprintf(STATUS, "%7.2f %7.2f %7.2f\n", native[a].xyz.x, native[a].xyz.y, native[a].xyz.z);
    //  fprintf(STATUS, "%5ld %5ld %5ld\n", native[a].xyz_int.x, native[a].xyz_int.y,
    //  native[a].xyz_int.z);

    if (data[a][b].check_clashes && distance < hard_core[native[a].smogtype][native[b].smogtype]) {
        data[a][b].clashes = data[b][a].clashes = 1;
        nclashes++;
    }

    if (data[a][b].check_contacts) {
        if ((distance <= contact_distance[native[a].smogtype][native[b].smogtype].b) &&
            (distance >= contact_distance[native[a].smogtype][native[b].smogtype].a)) {
            //	if(((native[a].smogtype==bb_O_type)&&(native[b].smogtype==bb_N_type))||
            //	   ((native[a].smogtype==bb_N_type)&&(native[b].smogtype==bb_O_type)))
            //	 {
            //          check_bb_contacts(a, b);
            //          fprintf(STATUS, "%4d %4d %4d %4d %2d %2d %2d %2d %2d\n",
            //	      a, b, native_residue[native[a].res_num].O,
            // native_residue[native[b].res_num].N, native[a].smogtype, native[b].smogtype,
            // native[a].res_num, native[b].res_num, abs(native[a].res_num-native[b].res_num));
            //	 }
            //        else
            //	 {
            data[a][b].contacts = data[b][a].contacts = 1;
            ncontacts++;
            //	 }
        }
        //    if((a==8)&&(b==497))
        //      fprintf(STATUS, "%5d %5d %10ld %10ld %10ld %5d %5d %5d\n",
    }
    return;
}

/* returns 0 if no changes */
/* returns 1 if only contacts changes */
/* returns 2 if only clashes changes */
/* returns 3 if both changes */

unsigned char CheckForDeltaContacts(struct contact_data *Data, struct int_vector XX,
                                    struct int_vector YY, short type_a, short type_b) {
    long int X_int, Y_int, Z_int;
    double   distance;

    X_int = XX.x - YY.x;
    Y_int = XX.y - YY.y;
    Z_int = XX.z - YY.z;

    distance = X_int * X_int + Y_int * Y_int + Z_int * Z_int;

    if (Data->check_contacts && Data->check_clashes) {
        if (((distance <= contact_distance[type_a][type_b].b) &&
             (distance >= contact_distance[type_a][type_b].a)) != Data->contacts) {
            Data->delta_contacts = !(Data->contacts);
            if ((distance < hard_core[type_a][type_b]) != (Data->clashes))
                Data->delta_clashes = mc_flags.clashed = !(Data->clashes);
            return 3;
        } else if ((distance < hard_core[type_a][type_b]) != (Data->clashes)) {
            Data->delta_clashes = mc_flags.clashed = !(Data->clashes);
            return 2;
        } else
            return 0;
    } else if (Data->check_clashes) {
        if ((distance < hard_core[type_a][type_b]) != Data->clashes) {
            Data->delta_clashes = mc_flags.clashed = !(Data->clashes);
            return 2;
        } else
            return 0;
    } else if (Data->check_contacts) {
        if (((distance <= contact_distance[type_a][type_b].b) &&
             (distance >= contact_distance[type_a][type_b].a)) != (Data->contacts)) {
            Data->delta_contacts = !(Data->contacts);
            return 1;
        } else
            return 0;
    }

    return 0;
}

//========================================================================================================
//========================================================================================================
void check_bb_contacts(short a, short b) {
    double        coo[3], oxy[3], cao[3];          // aceptor
    double        nit[3], con[3], can[3], hyd[3];  // donor
    double        len_ho, ang_h, ang_o, dih_coo;
    struct vector V3, H;
    int           coo_at, oxy_at, cao_at;
    int           nit_at, con_at, can_at;
    int           passed = 1;

    if ((native[a].smogtype == bb_O_type) && (native[b].smogtype == bb_N_type))  // a == O, b == N
    {
        if (b == 0)
            return;
        coo_at = native_residue[native[a].res_num].C;
        oxy_at = native_residue[native[a].res_num].O;
        cao_at = native_residue[native[a].res_num].CA;
        nit_at = native_residue[native[b].res_num].N;
        con_at = native_residue[native[b - 1].res_num].C;
        can_at = native_residue[native[b].res_num].CA;
    }

    else if ((native[a].smogtype == bb_N_type) &&
             (native[b].smogtype == bb_O_type))  // a ==O, b == N
    {
        if (a == 0)
            return;
        coo_at = native_residue[native[b].res_num].C;
        oxy_at = native_residue[native[b].res_num].O;
        cao_at = native_residue[native[b].res_num].CA;
        nit_at = native_residue[native[a].res_num].N;
        con_at = native_residue[native[a - 1].res_num].C;
        can_at = native_residue[native[a].res_num].CA;
    } else {
        fprintf(STATUS, "Error in bb_interactions!\n");
        exit(1);
    }

    coo[0] = native[coo_at].xyz.x;
    coo[1] = native[coo_at].xyz.y;
    coo[2] = native[coo_at].xyz.z;
    oxy[0] = native[oxy_at].xyz.x;
    oxy[1] = native[oxy_at].xyz.y;
    oxy[2] = native[oxy_at].xyz.z;
    nit[0] = native[nit_at].xyz.x;
    cao[0] = native[cao_at].xyz.x;
    cao[1] = native[cao_at].xyz.y;
    cao[2] = native[cao_at].xyz.z;
    nit[1] = native[nit_at].xyz.y;
    nit[2] = native[nit_at].xyz.z;
    con[0] = native[con_at].xyz.x;
    con[1] = native[con_at].xyz.y;
    con[2] = native[con_at].xyz.z;
    can[0] = native[can_at].xyz.x;
    can[1] = native[can_at].xyz.y;
    can[2] = native[can_at].xyz.z;
    MakeVector(native[nit_at].xyz, native[can_at].xyz, &H);
    MakeVector(native[nit_at].xyz, native[con_at].xyz, &V3);
    Add(V3, &H);
    Normalize(&H);
    Inverse(&H);
    Add(native[nit_at].xyz, &H);
    hyd[0] = H.x;
    hyd[1] = H.y;
    hyd[2] = H.z;

    c_bnd_len(hyd, oxy, &len_ho);
    c_bnd_ang(oxy, hyd, nit, &ang_h);
    c_bnd_ang(coo, oxy, hyd, &ang_o);
    c_dih_ang(cao, coo, oxy, hyd, &dih_coo);
    ang_h   = ang_h * rad2deg;
    ang_o   = ang_o * rad2deg;
    dih_coo = dih_coo * rad2deg;

    fprintf(STATUS, "%4d %4d %4d %4d %2d %2d %2d %2d %2d %7.3f %7.3f %7.3f %7.3f\n", a, b, oxy_at,
            nit_at, native[a].smogtype, native[b].smogtype, native[a].res_num, native[b].res_num,
            abs(native[a].res_num - native[b].res_num), len_ho, ang_h, ang_o, dih_coo);

    if (passed == 1) {
        data[a][b].contacts = data[b][a].contacts = 1;
        ncontacts++;
    }

    return;
}

void NewDeltaContacts(short rotate_natoms, short *rotate_atom, char *not_rotated) {
    int    i, j, k;
    short  temp;
    int    N, O, M, total_pairs = 0;
    short *A = 0;

    for (i = 0; i < rotate_natoms; i++) {
        N               = rotate_atom[i];
        temp_atom       = &native[N];
        temp            = temp_atom->smogtype;
        temp_xyz_int    = temp_atom->xyz_int;
        temp_cell3      = prev_native[N].matrix;
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
                        (data[N][M].check_contacts || data[N][M].check_clashes)) {
                        if (CheckForDeltaContacts(&data[N][M], temp_xyz_int, native[M].xyz_int,
                                                  temp, native[M].smogtype)) {
                            /* fprintf(STATUS, "PAIRS\t%d %d %d: %f\n", total_pairs, N, M,
                             * sqrt(D2(native[N].xyz, native[M].xyz)));*/
                            ab[total_pairs].a   = N;
                            ab[total_pairs++].b = M;
                        }
                        if (mc_flags.clashed)
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
                            if (is_rotated[N] != is_rotated[M] && data[N][M].check_contacts &&
                                data[N][M].contacts) {
                                /* ... then it is no longer in contact! */
                                /* so add pair to list, and default delta_contact = 0 turns off the
                                 * contact */
                                ab[total_pairs].a   = N;
                                ab[total_pairs++].b = M;
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
void fill_calpha_contact_string(const struct residue *residues, const struct atom *atoms,
                                short *contactstring) {
    int    i, j, atomi, atomj;
    double xi, yi, zi, xj, yj, zj, dist;
    for (i = 0; i < nresidues; i++) {
        atomi = residues[i].CA;
        xi    = atoms[atomi].xyz.x;
        yi    = atoms[atomi].xyz.y;
        zi    = atoms[atomi].xyz.z;
        for (j = i + min_seq_sep; j < nresidues; j++) /* CONTACTS DEF*/
        {
            atomj = residues[j].CA;
            xj    = atoms[atomj].xyz.x;
            yj    = atoms[atomj].xyz.y;
            zj    = atoms[atomj].xyz.z;
            dist  = ((xi - xj) * (xi - xj) + (yi - yj) * (yi - yj) + (zi - zj) * (zi - zj));
            if (dist < contact_calpha_cutoff * contact_calpha_cutoff) {
                contactstring[i * nresidues + j] = contactstring[j * nresidues + i] = 1;
            } else {
                contactstring[i * nresidues + j] = contactstring[j * nresidues + i] = 0;
            }
        }
    }
}

double number_of_calpha_contacts(const struct residue *residues, const struct atom *atoms) {
    int    i, j, atomi, atomj;
    double xi, yi, zi, xj, yj, zj, dist;
    double number_of_contacts = 0.;
    for (i = 0; i < nresidues; i++) {
        atomi = residues[i].CA;
        xi    = atoms[atomi].xyz.x;
        yi    = atoms[atomi].xyz.y;
        zi    = atoms[atomi].xyz.z;
        for (j = i + min_seq_sep; j < nresidues; j++) /* CONTACTS DEF*/
        {
            atomj = residues[j].CA;
            xj    = atoms[atomj].xyz.x;
            yj    = atoms[atomj].xyz.y;
            zj    = atoms[atomj].xyz.z;
            dist  = ((xi - xj) * (xi - xj) + (yi - yj) * (yi - yj) + (zi - zj) * (zi - zj));
            if (dist < contact_calpha_cutoff * contact_calpha_cutoff) {
                number_of_contacts += 1.;
            }
        }
    }
    return number_of_contacts;
}

double number_of_calpha_native_contacts(const struct residue *residues, const struct atom *atoms,
                                        const short *contactstring) {
    int    i, j, atomi, atomj;
    double xi, yi, zi, xj, yj, zj, dist;
    double number_of_contacts = 0.;
    for (i = 0; i < nresidues; i++) {
        atomi = residues[i].CA;
        xi    = atoms[atomi].xyz.x;
        yi    = atoms[atomi].xyz.y;
        zi    = atoms[atomi].xyz.z;
        for (j = i + min_seq_sep; j < nresidues; j++) /* CONTACTS DEF*/
        {
            atomj = residues[j].CA;
            xj    = atoms[atomj].xyz.x;
            yj    = atoms[atomj].xyz.y;
            zj    = atoms[atomj].xyz.z;
            dist  = ((xi - xj) * (xi - xj) + (yi - yj) * (yi - yj) + (zi - zj) * (zi - zj));
            if (dist < contact_calpha_cutoff * contact_calpha_cutoff &&
                contactstring[nresidues * i + j]) {
                number_of_contacts += 1.;
            }
        }
    }
    return number_of_contacts;
}

double number_of_calpha_nonnative_contacts(const struct residue *residues, const struct atom *atoms,
                                           const short *contactstring) {
    int    i, j, atomi, atomj;
    double xi, yi, zi, xj, yj, zj, dist;
    double number_of_contacts = 0.;
    for (i = 0; i < nresidues; i++) {
        atomi = residues[i].CA;
        xi    = atoms[atomi].xyz.x;
        yi    = atoms[atomi].xyz.y;
        zi    = atoms[atomi].xyz.z;
        for (j = i + min_seq_sep; j < nresidues; j++) /* CONTACTS DEF*/
        {
            atomj = residues[j].CA;
            xj    = atoms[atomj].xyz.x;
            yj    = atoms[atomj].xyz.y;
            zj    = atoms[atomj].xyz.z;
            dist  = ((xi - xj) * (xi - xj) + (yi - yj) * (yi - yj) + (zi - zj) * (zi - zj));
            if (dist < contact_calpha_cutoff * contact_calpha_cutoff &&
                !contactstring[nresidues * i + j]) {
                number_of_contacts += 1.;
            }
        }
    }
    return number_of_contacts;
}

double hamming_distance_calpha_contacts(const struct residue *residues, const struct atom *atoms,
                                        const short *contactstring) {
    int    i, j, atomi, atomj;
    double xi, yi, zi, xj, yj, zj, dist;
    short  iscontact;
    double hamming_distance = 0.;
    for (i = 0; i < nresidues; i++) {
        atomi = residues[i].CA;
        xi    = atoms[atomi].xyz.x;
        yi    = atoms[atomi].xyz.y;
        zi    = atoms[atomi].xyz.z;
        for (j = i + min_seq_sep; j < nresidues; j++) /* CONTACTS DEF*/
        {
            atomj     = residues[j].CA;
            xj        = atoms[atomj].xyz.x;
            yj        = atoms[atomj].xyz.y;
            zj        = atoms[atomj].xyz.z;
            dist      = ((xi - xj) * (xi - xj) + (yi - yj) * (yi - yj) + (zi - zj) * (zi - zj));
            iscontact = dist < contact_calpha_cutoff * contact_calpha_cutoff;
            hamming_distance += iscontact ^ contactstring[nresidues * i + j];
        }
    }
    return hamming_distance;
}
/* end wmj */