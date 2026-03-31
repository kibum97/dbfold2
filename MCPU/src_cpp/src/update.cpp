#include "update.h"
#include "init.h"
#include "lattice_util.h"
#include "vector.h"

void Update(
    struct Context *ctx,
    struct System  *sys,
    struct Simulation *sim
) { /* This function is called AFTER a move has been made. Counterintuitively, the
values for prev_E are set as the CURRETNT values of the energy. But that's because we make another
move before this function is called again, and so after that other move is made, the prev_E valious
really are previous

On the first step of the simulation, the function Reset_energies is used to compute the initial
energy, and the values for prev_E are set to that initial energy
*/
    int   i, j, k, temp1;
    short M, N, temp;

    for (i = 0; i < ctx->all_rotated_natoms; i++) {
        N              = ctx->all_rotated_atoms[i];
        temp_atom      = &ctx->native[N];
        temp_prev_atom = &ctx->prev_native[N];
        // identify matrix index j of the rotated atom N in the previous step
        j = 0;
        while (N != temp_prev_atom->matrix->atom_list[j])
            j++;
        // remove the rotated atom N from the previous matrix
        for (k = j; k < (temp_prev_atom->matrix->natoms - 1); k++)
            temp_prev_atom->matrix->atom_list[k] = temp_prev_atom->matrix->atom_list[k + 1];
        temp_prev_atom->matrix->natoms--;

        if (temp_atom->matrix->natoms) {
            j = 0;
            while (j < temp_atom->matrix->natoms && N > temp_atom->matrix->atom_list[j])
                j++;
            if (j != temp_atom->matrix->natoms) {
                temp                            = temp_atom->matrix->atom_list[j];
                temp_atom->matrix->atom_list[j] = N;
                for (k = j + 1; k < (temp_atom->matrix->natoms); k++) {
                    temp1                           = temp_atom->matrix->atom_list[k];
                    temp_atom->matrix->atom_list[k] = temp;
                    temp                            = temp1;
                }
                temp_atom->matrix->atom_list[temp_atom->matrix->natoms++] = temp;
            } else
                temp_atom->matrix->atom_list[temp_atom->matrix->natoms++] = N;

            if (temp_atom->matrix->natoms >=
                MAX_CELL_ATOMS) {  // We have some sort of problem if the number of atoms in the
                                   // matrix's atom_list exceeds or is equal to MAX_CELL_ATOMS (100)
                fprintf(
                    sim->STATUS,
                    "Lattice Error: Update(), num. of atoms exceeds maximum num. %d of the cell",
                    MAX_CELL_ATOMS);
                fprintf(sim->STATUS, "Lattice Error: atom %4s %4d %4s\n", temp_atom->atomname,
                        temp_atom->res_num, temp_atom->res);
                fprintf(sim->STATUS, "Lattice Error: atom lists\n");
                for (k = 0; k < temp_atom->matrix->natoms;
                     k++) {  // Loop through all entries in the atom_list for the matrix affiliated
                             // with the current atom
                    M = temp_atom->matrix
                            ->atom_list[k];  // The value for kth entry in atom list, indexes some
                                             // other atom in the native structure
                    fprintf(sim->STATUS, "Lattice Error: %4d %4s %4d %4s\n", M, ctx->native[M].atomname,
                            ctx->native[M].res_num, ctx->native[M].res);
                }
                exit(1);
            }
        } else
            temp_atom->matrix->atom_list[temp_atom->matrix->natoms++] = N;

        Copy(temp_atom->xyz, &temp_prev_atom->xyz);
        CopyLatticeCoordinates(ctx->native[N], &ctx->prev_native[N]);
    }

    ctx->E_pot += ctx->dE_pot;
    ctx->E_tor += ctx->dE_tor;
    ctx->E_sct += ctx->dE_sct;
    ctx->E_aro += ctx->dE_aro;
    ctx->E_hbond += ctx->dE_hbond;
    ctx->E_constraint += ctx->dE_constraint;
    ctx->E += ctx->dE;
    ctx->prev_E            = ctx->E;
    ctx->prev_E_pot        = ctx->E_pot;
    ctx->prev_E_tor        = ctx->E_tor;
    ctx->prev_E_sct        = ctx->E_sct;
    ctx->prev_E_aro        = ctx->E_aro;
    ctx->prev_E_hbond      = ctx->E_hbond;
    ctx->prev_E_constraint = ctx->E_constraint;
    // fprintf(stdout, "Update : dE_sct : %.5f, E :%.5f, dE : %.5f, prev_E : %.5f\n", dE_sct, E, dE,
    // prev_E); fprintf(stdout, "%9.5f %9.5f\n", E, dE);

    ctx->ncontacts += ctx->delta_contacts;

    //  if ((sidechain_step!=0)||(USE_ROTAMERS))
    if (ctx->sidechain_step != 0)
        for (i = 0; i < ctx->native_residue[ctx->mc.sel_res_num].ntorsions; i++)
            ctx->native_residue[ctx->mc.sel_res_num].chi[i] += ctx->mc.delta_angle[i];
    else
        for (i = 0; i < ctx->mc.loop_size; i++) {
            ctx->native_residue[ctx->mc.selected[i]].phi += ctx->mc.delta_phi_angle[i];
            ctx->native_residue[ctx->mc.selected[i]].psi += ctx->mc.delta_psi_angle[i];
        }

    for (i = 0; i < ctx->total_pairs; i++) {
        M                   = ctx->ab[i].a;
        N                   = ctx->ab[i].b;
        ctx->data[M][N].contacts = ctx->data[N][M].contacts = ctx->data[M][N].delta_contacts;
        ctx->data[M][N].delta_contacts                 = 0;
        if (!ctx->mc_flags.init)
            ctx->data[M][N].delta_clashes = 0;
    }

    if (!ctx->mc_flags.init) {
        for (i = 0; i < ctx->total_hbond_pairs; i++) {
            M                = ctx->hbond_pair[i].a;
            N                = ctx->hbond_pair[i].b;
            ctx->data[M][N].hbond = ctx->data[N][M].hbond = ctx->data[M][N].delta_hbond;
            ctx->data[M][N].closehb = ctx->data[N][M].closehb = ctx->data[M][N].delta_closehb;
            ctx->data[M][N].delta_hbond = ctx->data[N][M].delta_hbond = sys->NO_HBOND;
            ctx->data[M][N].delta_closehb = ctx->data[N][M].delta_closehb = 3;
        }
    }

    if (ctx->mc_flags.init) {
        ctx->nclashes += ctx->delta_nclashes;
        for (i = 0; i < ctx->total_pairs2; i++) {
            M                  = ctx->cd[i].a;
            N                  = ctx->cd[i].b;
            ctx->data[M][N].clashes = ctx->data[N][M].clashes = ctx->data[M][N].delta_clashes;
            ctx->data[M][N].delta_clashes                = 0;
        }
    }
    return;
}

void Restore(
    struct Context *ctx,
    struct MCIntegrator *integrator,
    struct System  *sys
) {
    int   i;
    short M, N;

    // fprintf(stdout, "Restore(): mc.sel_res_num %d, sidechain_step %1d, sidemovedone %1d\n",
    // mc.sel_res_num, sidechain_step, sidemovedone);

    if ((ctx->sidechain_step != 0) && (ctx->sidemovedone != 0))
        for (i = 0; i < ctx->native_residue[ctx->mc.sel_res_num].ntorsions; i++)
            ctx->native_residue[ctx->mc.sel_res_num].tmpchi[i] -= ctx->mc.delta_angle[i];

    for (i = 0; i < ctx->all_rotated_natoms; i++) {
        N              = ctx->all_rotated_atoms[i];
        temp_atom      = &ctx->native[N];
        temp_prev_atom = &ctx->prev_native[N];
        Copy(temp_prev_atom->xyz, &temp_atom->xyz);
        CopyLatticeCoordinates(ctx->prev_native[N], &ctx->native[N]);
    }

    if (integrator->USE_ROTAMERS)
        ctx->cur_rotamers[ctx->mc.sel_res_num] = ctx->old_rotamer;

    for (i = 0; i < ctx->total_pairs; i++) {
        M                         = ctx->ab[i].a;
        N                         = ctx->ab[i].b;
        ctx->data[M][N].delta_contacts = 0;
        ctx->data[M][N].delta_clashes  = 0;
    }

    if (ctx->mc_flags.init) {
        for (i = 0; i < ctx->total_pairs2; i++) {
            M                        = ctx->cd[i].a;
            N                        = ctx->cd[i].b;
            ctx->data[M][N].delta_clashes = 0;
        }
    }

    else {
        for (i = 0; i < ctx->total_hbond_pairs; i++) {
            M                      = ctx->hbond_pair[i].a;
            N                      = ctx->hbond_pair[i].b;
            ctx->data[M][N].delta_hbond = ctx->data[N][M].delta_hbond = sys->NO_HBOND;
            ctx->data[M][N].delta_closehb = ctx->data[N][M].delta_closehb = 3;
        }
    }

    ctx->E            = ctx->prev_E;
    ctx->E_pot        = ctx->prev_E_pot;
    ctx->E_tor        = ctx->prev_E_tor;
    ctx->E_sct        = ctx->prev_E_sct;
    ctx->E_aro        = ctx->prev_E_aro;
    ctx->E_hbond      = ctx->prev_E_hbond;
    ctx->E_constraint = ctx->prev_E_constraint;

    return;
}
