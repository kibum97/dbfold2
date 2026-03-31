#include "move.h"  // Already there

#include "constraint.h"
#include "contacts.h"  // For number_of_calpha_native_contacts
#include "define.h"    // For constants like CLUSTER_MOVE
#include "energy.h"    // For torsionenergy(), sctenergy(), aromaticenergy()
#include "getrms.h"
#include "globals.h"  // For natives, nresidues, constraint_align, etc.
#include "hbonds.h"
#include "init.h"
#include "lattice_util.h"  // For UpdateLattice(), Update(), Restore()
#include "loop.h"          // For integloop()
#include "misc_util.h"
#include "rng.h"  // For threefryrand(), GaussianNum()
#include "rotate.h"
#include "update.h"
#include "vector.h"  // For DoRotation()

// You also have errors regarding string comparisons and standard I/O
#include <stdio.h>   // For potential fprintf/sprintf calls
#include <string.h>  // For strcmp

/*==================================================*/
/*                  main program                    */
/*==================================================*/

double new_rms;

void MakeMove(struct Context *ctx, struct System *sys, struct Topology *top,
              struct Simulation *sim, struct MCIntegrator *integrator, float step_size,
              float use_global_bb_moves) {
    int kkk;
    int reject, del, N, M, i;
    ctx->sidechain_step = 0;
    int use_yang        = 0;
    int n_soln          = 0;
    //  float hb_before_move = 0.;
    //  float hb_after_move = 0.;
    //  int j;
    //  float e;
    struct alignment align;

    align.NFRAG           = 1;
    align.seqptr[1].x1    = 1;
    align.structptr[1].x1 = 1;
    align.seqptr[1].x2    = top->nresidues;
    align.structptr[1].x2 = top->nresidues;

    struct alignment constraint_align = align;

    ctx->backbone_accepted = 0;

    /* wmj ************************************* */
    ctx->diff_number_of_contacts_current = ctx->natives - sim->number_of_contacts_setpoint;

    /* end wmj ********************************* */

    do { /* See the last line of this clause, while (sidechain_step++ < SIDECHAIN_MOVES)
          In other words, you make a backbone move first (since sidechain_step above is set to 0
          Then you make however many sidechain moves as dictated by variable SIDECHAIN_MOVES
          i.e. with SIDECHAIN_MOVES = 1, the loop execs twice...One backbone, one sidechain  */

        reject                = 0;
        ctx->mc_flags.clashed = 0;
        sim->nomove           = 0;
        ctx->sidemovedone     = 0;

        /* make backbone move */
        if (ctx->sidechain_step == 0) {
            if ((use_global_bb_moves) && (use_global_bb_moves > 0.001)) {
                fprintf(sim->STATUS, "USE_GLOBAL_BB_MOVES is turned on!!!");
                exit(1);
            }

            // Try one of the various possible backbone moves
            if ((threefryrand() < integrator->CLUSTER_MOVE) &&
                integrator->USE_CLUSTER) {  // Cluster_move will often be set to 0 to suppress
                                            // knowledge moves
                // fprintf(STATUS, "LoopBackbondMove():\n");
                LoopBackboneMove(ctx, sim, integrator, top, sys,
                                 step_size);  // Make a knowledge based backbone move
                use_yang = 0;
            } else {
                if (integrator->YANG_MOVE &&
                    threefryrand() <
                        integrator->YANG_MOVE) {    // Make a local move, as in Dill paper
                    integloop(ctx, top, integrator, sys, sim, step_size, &n_soln);  // See loop.h
                    use_yang = 1;
                } else {
                    // fprintf(STATUS, "LocalBackboneMove():\n");
                    LocalBackboneMove(
                        ctx, top, sys, sim, integrator,
                        step_size);  // Counterintuitively, this is actually a global move
                    exit(1);
                    use_yang = 0;
                }
            } /* end (sidechain_step != 0) */

            /* VZ
             * struct_f2, the structure backbone, is updated here
             * because a backbone move was attempted.
             * Acts only if rmsd_constraint is active
             * Note: struct_f1 and struct_f2 are declared in backbone.h,
             *   and struct_f1 is set in fold.h
             */
            if (sys->rmsd_constraint > 0) {
                for (i = 0; i < top->nalign; ++i) {
                    ctx->struct_f2[i + 1].CA.x =
                        ctx->native[ctx->native_residue[top->map_to_seq[i]].CA].xyz.x;
                    ctx->struct_f2[i + 1].CA.y =
                        ctx->native[ctx->native_residue[top->map_to_seq[i]].CA].xyz.y;
                    ctx->struct_f2[i + 1].CA.z =
                        ctx->native[ctx->native_residue[top->map_to_seq[i]].CA].xyz.z;
                }
                new_rms = getrms(ctx->struct_f1, ctx->struct_f2, constraint_align);
            }
        } else { /* make sidechain move (sidechain_step != 0) */
            if (top->total_ntorsions != 0) {
                // fprintf(STATUS, "SidechainMove():\n");
                SidechainMove(ctx, top, integrator, sys);
            }
            use_yang = 0;
        }

        if (!ctx->mc_flags.init && ctx->mc_flags.clashed) {
            reject = 1;
            if ((sim->nomove == 0) && (ctx->sidechain_step == 0))
                sim->nrejected++;
        } else if ((use_yang == 1) && (n_soln == 0)) {
            reject = 1;
            // fprintf(STATUS, "Yang move rejected \n");
            if (ctx->sidechain_step == 0)
                sim->nothers++;
        }
        // else if ((use_yang == 1) && (n_soln != 0)) {
        // fprintf(STATUS, "Yang move accepted \n");

        //}
        /* below block is new code to support RMSD constraint
         * commenting this block fixes bug.
         * --> has to be something messing with logic?
         */
        else if ((ctx->sidechain_step == 0) && (sys->rmsd_constraint > 0) &&
                 (new_rms > sys->rmsd_constraint)) {
            reject = 1;
            sim->nrejected++;
        } else {
            ctx->delta_nclashes = 0;
            if (ctx->mc_flags.init)
                for (i = 0; i < ctx->total_pairs2; i++)
                    ctx->delta_nclashes += ctx->data[ctx->cd[i].a][ctx->cd[i].b].delta_clashes -
                                           ctx->data[ctx->cd[i].a][ctx->cd[i].b].clashes;

            ctx->dE            = 0;
            ctx->dE_pot        = 0;
            ctx->dE_hbond      = 0;
            ctx->dE_tor        = 0;
            ctx->dE_aro        = 0;
            ctx->dE_sct        = 0;
            ctx->dE_constraint = 0;

            if (ctx->sidechain_step == 0)
                ctx->dE_tor = torsionenergy(ctx, top, sys) - ctx->prev_E_tor;

            if ((ctx->sidechain_step != 0) && ctx->sidemovedone)
                ctx->dE_sct = sctenergy(ctx, sys, top) - ctx->prev_E_sct;

            ctx->dE_aro         = aromaticenergy(ctx, sys, top) - ctx->prev_E_aro;
            ctx->delta_contacts = 0;
            for (i = 0; i < ctx->total_pairs; i++) { /*AB: Loop through all pairs of atoms and
                                                   compute potential change due to move*/
                N   = ctx->ab[i].a;
                M   = ctx->ab[i].b; /*N, M are indices for pairs of atoms*/
                del = ctx->data[N][M].delta_contacts -
                      ctx->data[N][M].contacts; /*I think del indicates somethign like how many new
                                              contacts have formed?*/
                ctx->delta_contacts += del;
                ctx->dE_pot +=
                    ((float)del) *
                    sys->potential[ctx->native[N].smogtype]
                                  [ctx->native[M].smogtype]; /*Compute change in mu potential
                                                           contribution ...I think the del in front
                                                           somehow restricts your attention to
                                                           atoms that have come into contact?*/
            }
            if (sys->weight_hbond) {
                if (ctx->sidechain_step == 0)
                    ctx->dE_hbond = FoldHydrogenBonds(top, ctx, sys) - ctx->prev_E_hbond;
            }

            /* here, code used to test for rmsd */
            //	align_drms(native, native_residue, struct_native, struct_residue, map_to_seq,
            // map_to_struct, nalign, &bb_rms);

            if (strcmp(sim->constraint_file, "None") != 0) {  // AB
                ctx->dE_constraint = Compute_constraint_energy(ctx, sys, ctx->native_residue, ctx->native) -
                                     ctx->prev_E_constraint;  // AB
            } else if (strcmp(sim->rmsd_constraint_file, "None") != 0) {
                ctx->dE_constraint = Compute_frag_rmsd_constraint_energy(
                                         sys, ctx->struct_f1, ctx->struct_f2, constraint_align) -
                                     ctx->prev_E_constraint;  // AB
            }  // AB: Otherwise it is 0, as previously initialized

            ctx->dE = sys->weight_potential * ctx->dE_pot +
                      sys->weight_clash * ctx->delta_nclashes + sys->weight_hbond * ctx->dE_hbond +
                      TOR_WEIGHT * ctx->dE_tor + ARO_WEIGHT * ctx->dE_aro +
                      SCT_WEIGHT * ctx->dE_sct + ctx->dE_constraint;  // AB addded last term
            double arg = ctx->dE / integrator->MC_TEMP;
            /* wmj ****************************** */

            ctx->new_natives = number_of_calpha_native_contacts(
                top, ctx->native_residue, ctx->native, ctx->orig_contactstring);
            ctx->diff_number_of_contacts_new = ctx->new_natives - sim->number_of_contacts_setpoint;
            double dbias =
                sys->k_bias *
                (ctx->diff_number_of_contacts_new * ctx->diff_number_of_contacts_new -
                 ctx->diff_number_of_contacts_current * ctx->diff_number_of_contacts_current);
            arg += dbias;
            /* end wmj ************************** */

            /************************  Start AB ***********************************/
            double acceptance_crit = exp(-arg);

            // Jacobi fix implemented 7/17 by AB
            if (use_yang == 1) {
                acceptance_crit = acceptance_crit * (n_soln / (double)ctx->soln_no_before) *
                                  (ctx->jacobi_after /
                                   ctx->jacobi_before);  // AB changed to obey detailed balance in
                                                         // case local moves are used
                // fprintf(STATUS, "jacobi_after = %7.3f, jacobi_before = %7.3f, n_solutions_after =
                // %i, n_solutions_before = %i, dE/MC_TEMP = %7.3f, acceptance criterion is %7.5f
                // \n", jacobi_after, jacobi_before, n_soln, soln_no_before,arg, acceptance_crit);
                // if (acceptance_crit>=1 || threefryrand() < acceptance_crit){
                // fprintf(STATUS, "Move may have been accepted \n");
                // };
            }
            /*End AB*/

            /* Metropolis */

            /*AB commented out: */
            // if (arg <= 0 || threefryrand() < exp(-arg)){

            /*AB replaced the above with*/
            if (acceptance_crit >= 1 || threefryrand() < acceptance_crit) {  // Move is accepted!
                // if (use_yang ==1){
                // fprintf(STATUS, "Move accepted \n");
                // }

                // end AB
                Update(ctx, sys, sim);  // Update() is where you actually update the coordinates of the structure
                for (kkk = 0; kkk < sys->N_constraints; kkk++) {
                    ctx->disulfide_pairs[kkk] = ctx->disulfide_pairs_attempt[kkk];
                }
                ctx->natives = ctx->new_natives;
                // fprintf(STATUS, "Updated ...\n");
                // ResetEnergies();
                if ((sim->nomove == 0) && (ctx->sidechain_step == 0)) {
                    sim->naccepted++;
                    ctx->backbone_accepted = 1;
                } else
                    sim->n_sidechain_accepted++;
                //    fprintf(STATUS, "%7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n", E,
                //    0.2*E_pot+0.8*E_hbond, E_pot, E_hbond, dE, dE_pot, dE_hbond);
            } else {
                reject = 1;
                if ((sim->nomove == 0) && (ctx->sidechain_step == 0))
                    sim->nrejected++;
            }
        }

        if (reject) {
            Restore(ctx, integrator, sys);
            // fprintf(STATUS, "rejected ...\n");
            // ResetEnergies();
        }

        for (i = 0; i < ctx->all_rotated_natoms; i++) {
            is_rotated[ctx->all_rotated_atoms[i]] = 0;
        }

    } while (ctx->sidechain_step++ < integrator->SIDECHAIN_MOVES);

    return;
}

/*==================================================*/

void LocalBackboneMove(struct Context *ctx, struct Topology *top, struct System *sys,
                       struct Simulation *sim, struct MCIntegrator *integrator, float step_size) {
    ctx->mc.loop_size       = 1;
    ctx->all_rotated_natoms = 0;
    ctx->total_pairs = ctx->total_pairs2 = 0;
    ctx->total_hbond_pairs               = 0;
    sim->nomove                          = 0;

    ctx->mc.sel_res_num = (int)(threefryrand() * top->nresidues);
    if (top->is_template[ctx->mc.sel_res_num] == 1) {
        sim->nomove = 1;
        sim->nothers++;
        return;
    }
    if (ctx->native_residue[ctx->mc.sel_res_num].amino_num ==
        14) { /* cannot rotate around proline phi */
        if (ctx->mc.sel_res_num != top->nresidues - 1)
            ctx->mc.is_phi = 0;
        else { /* if proline is last residue, no rotation possible */
            sim->nomove = 1;
            sim->nothers++;
            return;
        }
    } else if (ctx->mc.sel_res_num == 0)
        ctx->mc.is_phi = 0;
    else if (ctx->mc.sel_res_num == top->nresidues - 1)
        ctx->mc.is_phi = 1;
    else
        ctx->mc.is_phi = (int)(threefryrand() * 2);
    // fprintf(STATUS, "mainchain move at %5d\n", ctx->mc.sel_res_num);

    step_size = step_size * GaussianNum();
    ;

    BackboneMove(ctx, top, integrator, step_size);
    UpdateLattice(ctx, integrator, integrator->rotate_natoms[ctx->mc.is_phi][ctx->mc.sel_res_num],
                  integrator->rotate_atom[ctx->mc.is_phi][ctx->mc.sel_res_num]);

    for (int i = 0; i < integrator->rotate_natoms[ctx->mc.is_phi][ctx->mc.sel_res_num]; i++) {
        is_rotated[integrator->rotate_atom[ctx->mc.is_phi][ctx->mc.sel_res_num][i]] = 1;
    }

    NewDeltaContacts(ctx, sys, integrator->rotate_natoms[ctx->mc.is_phi][ctx->mc.sel_res_num],
                     integrator->rotate_atom[ctx->mc.is_phi][ctx->mc.sel_res_num],
                     integrator->not_rotated[ctx->mc.is_phi][ctx->mc.sel_res_num]);
    return;
}

void LoopBackboneMove(struct Context *ctx, struct Simulation *sim,
                      struct MCIntegrator *integrator, struct Topology *top,
                      struct System *sys, float absolute_step_size) {
    /*Contrary to the naming, this one actually does the knowledge based move*/

    int   a, b, c, d, i, j;
    float step_phi, step_psi;
    float desire_phi, desire_psi;
    int   cluster_bin;
    float use_cluster;

    /* every 1 moves, a new triplet of residues is selected. */
    /* each move, the phi/psi angles of a single residue are changed */

    ctx->mc.loop_size = 1;

    /* triplet is selected */
    if (0 == 0) {
        ctx->mc.sel_triplet = (int)(threefryrand() * integrator->TOTAL_TRIPLE_LOOP_MOVES) +
                              integrator->TOTAL_SINGLE_LOOP_MOVES +
                              integrator->TOTAL_DOUBLE_LOOP_MOVES;

        if (integrator->residue_triplets[ctx->mc.sel_triplet].a > top->nresidues / 2.0) {
            ctx->mc.selected[0] = integrator->residue_triplets[ctx->mc.sel_triplet].a;
            ctx->mc.selected[1] = integrator->residue_triplets[ctx->mc.sel_triplet].b;
            ctx->mc.selected[2] = integrator->residue_triplets[ctx->mc.sel_triplet].c;
        } else {
            if (integrator->residue_triplets[ctx->mc.sel_triplet].c >= 0) {
                ctx->mc.selected[0] = integrator->residue_triplets[ctx->mc.sel_triplet].c;
                ctx->mc.selected[1] = integrator->residue_triplets[ctx->mc.sel_triplet].b;
                ctx->mc.selected[2] = integrator->residue_triplets[ctx->mc.sel_triplet].a;
            } else if (integrator->residue_triplets[ctx->mc.sel_triplet].b >= 0) {
                ctx->mc.selected[0] = integrator->residue_triplets[ctx->mc.sel_triplet].b;
                ctx->mc.selected[1] = integrator->residue_triplets[ctx->mc.sel_triplet].a;
                ctx->mc.selected[2] = integrator->residue_triplets[ctx->mc.sel_triplet].c;
            } else {
                ctx->mc.selected[0] = integrator->residue_triplets[ctx->mc.sel_triplet].a;
                ctx->mc.selected[1] = integrator->residue_triplets[ctx->mc.sel_triplet].b;
                ctx->mc.selected[2] = integrator->residue_triplets[ctx->mc.sel_triplet].c;
            }
        }
    }
    if (top->is_template[ctx->mc.selected[0]] == 1) {
        sim->nomove = 1;
        sim->nothers++;
        return;
    }

    ctx->all_rotated_natoms = 0;
    ctx->total_pairs = ctx->total_pairs2 = 0;
    ctx->total_hbond_pairs               = 0;

    /* move-cycle defines the residue to be changed */
    i               = ctx->move_cycle % 1;
    ctx->move_cycle = i;
    //  ++move_cycle;

    check_phipsi(top, ctx);
    cluster_bin = (int)(NOCLUSTERS * threefryrand());
    //  mc.selected[i] = 92;
    //  cluster_bin = 1;
    desire_phi = sys->cluster_phi[ctx->native_residue[ctx->mc.selected[i]].amino_num][cluster_bin];
    desire_psi = sys->cluster_psi[ctx->native_residue[ctx->mc.selected[i]].amino_num][cluster_bin];
    // fprintf(STATUS, "%d\n", cluster_bin);
    /* rotates phi and psi */
    if (top->secstr[ctx->mc.selected[i]] == 'H')
        use_cluster = 0.5;
    else if (top->secstr[ctx->mc.selected[i]] == 'E')
        use_cluster = 0.0;
    else if (top->secstr[ctx->mc.selected[i]] == 'L')
        use_cluster = 0.0;
    else if (top->secstr[ctx->mc.selected[i]] == 'C')
        use_cluster = integrator->USE_CLUSTER;
    else {
        fprintf(sim->STATUS, "ERROR! secondary structure prediction has a unknown charactrer: %c\n",
                top->secstr[ctx->mc.selected[i]]);
        exit(1);
    }
    if (ctx->mc.selected[i] > top->nresidues / 2.0) {
        if (ctx->mc.selected[i] != 0 && ctx->mc.selected[i] != (top->nresidues - 1) &&
            ctx->native_residue[ctx->mc.selected[i]].amino_num != 14) {
            /* can't rotate around a proline phi */
            if (use_cluster > threefryrand()) {
                step_phi = desire_phi - ctx->cur_phi[ctx->mc.selected[i] - 1];
                step_phi += GaussianNum() * CLUSTER_NOISE;
                // fprintf(STATUS, "Did a knowledge based move at step %10ld !!\n", mcstep);
            } else
                // fprintf(STATUS, "LoopMove called, but no knowledge based move done\n");
                step_phi = GaussianNum() * CLUSTER_NOISE;
            step_phi *= deg2rad;
            //	fprintf(STATUS, "%f\n", step_size);

            a = ctx->native_residue[ctx->mc.selected[i] - 1].C;
            b = ctx->native_residue[ctx->mc.selected[i]].N;
            c = ctx->native_residue[ctx->mc.selected[i]].CA;
            d = ctx->native_residue[ctx->mc.selected[i]].C;
            DoRotation(ctx, a, b, c, d, step_phi, integrator->loop_rotate_natoms[ctx->mc.selected[i]][0],
                       integrator->loop_rotate_atoms[ctx->mc.selected[i]][0]);
            ctx->mc.delta_phi_angle[i] = step_phi;
            for (j = 0; j < integrator->loop_rotate_natoms[ctx->mc.selected[i]][0]; j++) {
                is_rotated[integrator->loop_rotate_atoms[ctx->mc.selected[i]][0][j]] = 1;
            }
        }
        if (ctx->mc.selected[i] != 0 && ctx->mc.selected[i] != (top->nresidues - 1)) {
            if (use_cluster > threefryrand()) {
                step_psi = desire_psi - ctx->cur_psi[ctx->mc.selected[i] - 1];
                step_psi += GaussianNum() * CLUSTER_NOISE;
                // fprintf(STATUS, "Did a knowledge based move at step %10ld !!\n", mcstep);
                // fprintf(STATUS, "Did a knowledge based move!!\n", secstr[mc.selected[i]]);
            } else
                // fprintf(STATUS, "LoopMove called, but no knowledge based move done\n");
                step_psi = GaussianNum() * CLUSTER_NOISE;
            step_psi *= deg2rad;

            a = ctx->native_residue[ctx->mc.selected[i]].N;
            b = ctx->native_residue[ctx->mc.selected[i]].CA;
            c = ctx->native_residue[ctx->mc.selected[i]].C;
            d = ctx->native_residue[ctx->mc.selected[i] + 1].N;
            DoRotation(ctx, a, b, c, d, step_psi, integrator->loop_rotate_natoms[ctx->mc.selected[i]][1],
                       integrator->loop_rotate_atoms[ctx->mc.selected[i]][1]);
            ctx->mc.delta_psi_angle[i] = step_psi;
            for (j = 0; j < integrator->loop_rotate_natoms[ctx->mc.selected[i]][1]; j++) {
                is_rotated[integrator->loop_rotate_atoms[ctx->mc.selected[i]][1][j]] = 2;
            }
        }
    } else {
        if (ctx->mc.selected[i] != 0 && ctx->mc.selected[i] != (top->nresidues - 1)) {
            if (use_cluster > threefryrand()) {
                step_psi = desire_psi - ctx->cur_psi[ctx->mc.selected[i] - 1];
                step_psi += GaussianNum() * CLUSTER_NOISE;
                // fprintf(STATUS, "Did a knowledge based move at step %10ld !!\n", mcstep);
            } else
                step_psi = GaussianNum() * CLUSTER_NOISE;
            step_psi *= -deg2rad;

            a = ctx->native_residue[ctx->mc.selected[i]].N;
            c = ctx->native_residue[ctx->mc.selected[i]].CA;
            b = ctx->native_residue[ctx->mc.selected[i]].C;
            d = ctx->native_residue[ctx->mc.selected[i] + 1].N;
            DoRotation(ctx, a, b, c, d, -step_psi,
                       integrator->loop_rotate_natoms[ctx->mc.selected[i]][0],
                       integrator->loop_rotate_atoms[ctx->mc.selected[i]][0]);
            ctx->mc.delta_psi_angle[i] = -step_psi;
            for (j = 0; j < integrator->loop_rotate_natoms[ctx->mc.selected[i]][0]; j++) {
                is_rotated[integrator->loop_rotate_atoms[ctx->mc.selected[i]][0][j]] = 1;
            }
        }
        if (ctx->mc.selected[i] != 0 && ctx->mc.selected[i] != (top->nresidues - 1) &&
            ctx->native_residue[ctx->mc.selected[i]].amino_num != 14) {
            if (use_cluster > threefryrand()) {
                step_phi = desire_phi - ctx->cur_phi[ctx->mc.selected[i] - 1];
                step_phi += GaussianNum() * CLUSTER_NOISE;
            } else
                step_phi = GaussianNum() * CLUSTER_NOISE;
            step_phi *= -deg2rad;

            a = ctx->native_residue[ctx->mc.selected[i] - 1].C;
            c = ctx->native_residue[ctx->mc.selected[i]].N;
            b = ctx->native_residue[ctx->mc.selected[i]].CA;
            d = ctx->native_residue[ctx->mc.selected[i]].C;
            DoRotation(ctx, a, b, c, d, -step_phi,
                       integrator->loop_rotate_natoms[ctx->mc.selected[i]][1],
                       integrator->loop_rotate_atoms[ctx->mc.selected[i]][1]);
            ctx->mc.delta_phi_angle[i] = -step_phi;
            for (j = 0; j < integrator->loop_rotate_natoms[ctx->mc.selected[i]][1]; j++) {
                is_rotated[integrator->loop_rotate_atoms[ctx->mc.selected[i]][1][j]] = 2;
            }
        }
    }
    //  check_phipsi();
    //   fprintf(STATUS, "%3d %s %3d %3d %7.5f %7.5f %7.5f %7.5f\n",
    //       mc.selected[i], native_residue[mc.selected[i]].res,
    //       native_residue[mc.selected[i]].amino_num, cluster_bin, desire_phi,
    //       cur_phi[mc.selected[i]-1], desire_psi, cur_psi[mc.selected[i]-1]);

    UpdateLattice(ctx, integrator, integrator->loop_rotate_natoms[ctx->mc.selected[i]][0],
                  integrator->loop_rotate_atoms[ctx->mc.selected[i]][0]);
    ctx->all_rotated_natoms = integrator->loop_rotate_natoms[ctx->mc.selected[i]][0];
    ctx->all_rotated_atoms  = integrator->loop_rotate_atoms[ctx->mc.selected[i]][0];

    NewDeltaContacts(ctx, sys, integrator->loop_rotate_natoms[ctx->mc.selected[i]][0],
                     integrator->loop_rotate_atoms[ctx->mc.selected[i]][0],
                     integrator->loop_not_rotated[ctx->mc.selected[i]][0]);
    return;
}

void BackboneMove(struct Context *ctx, struct Topology *top, struct MCIntegrator *integrator,
                  float step_size) {
    short temp;
    int   a, b, c, d; /*          d    */
                      /*         /     */
                      /*    b - c      */
                      /*   /           */
                      /*  a            */

    ctx->mc.delta_angle[0] = step_size;

    if (ctx->mc.is_phi) {
        a                          = ctx->native_residue[ctx->mc.sel_res_num - 1].C;
        b                          = ctx->native_residue[ctx->mc.sel_res_num].N;
        c                          = ctx->native_residue[ctx->mc.sel_res_num].CA;
        d                          = ctx->native_residue[ctx->mc.sel_res_num].C;
        ctx->mc.delta_phi_angle[0] = step_size; /* just for proper updating purposes */
        ctx->mc.delta_psi_angle[0] = 0;
    } else {
        a                          = ctx->native_residue[ctx->mc.sel_res_num].N;
        b                          = ctx->native_residue[ctx->mc.sel_res_num].CA;
        c                          = ctx->native_residue[ctx->mc.sel_res_num].C;
        d                          = ctx->native_residue[ctx->mc.sel_res_num + 1].N;
        ctx->mc.delta_psi_angle[0] = step_size;
        ctx->mc.delta_phi_angle[0] = 0;
    }

    if (ctx->mc.sel_res_num <= top->nresidues / 2.0) {
        ctx->mc.delta_angle[0] *= -1;
        ctx->mc.delta_psi_angle[0] *= -1;
        ctx->mc.delta_phi_angle[0] *= -1;
        temp = b;
        b    = c;
        c    = temp;
    }

    DoRotation(ctx, a, b, c, d, ctx->mc.delta_angle[0],
               integrator->rotate_natoms[ctx->mc.is_phi][ctx->mc.sel_res_num],
               integrator->rotate_atom[ctx->mc.is_phi][ctx->mc.sel_res_num]);

    ctx->all_rotated_natoms = integrator->rotate_natoms[ctx->mc.is_phi][ctx->mc.sel_res_num];
    ctx->all_rotated_atoms  = integrator->rotate_atom[ctx->mc.is_phi][ctx->mc.sel_res_num];

    return;
}

void MakeSidechainMove(struct Context *ctx, struct Topology *top,
                       struct MCIntegrator *integrator, struct System *sys) {
    int   a, b, c, d, i, j;
    int   sel_no_chi;
    float p_0to1;
    float cummul_prob;

    ctx->mc.sel_rotamer =
        (int)(threefryrand() * ctx->native_residue[ctx->mc.sel_res_num].nrotamers);
    if (integrator->USE_ROT_PROB == 1) {
        p_0to1      = threefryrand() * 100;
        cummul_prob = 0.;
        for (i = 0; i < top->no_chi_list[ctx->native_residue[ctx->mc.sel_res_num].amino_num]; i++) {
            cummul_prob += sys->prob_ang[ctx->native_residue[ctx->mc.sel_res_num].amino_num][i];
            if (cummul_prob > p_0to1)
                break;
        }
        sel_no_chi = i;
    } else
        sel_no_chi = (int)(threefryrand() *
                           top->no_chi_list[ctx->native_residue[ctx->mc.sel_res_num].amino_num]);
    // fprintf(STATUS, "sidechain move at %5d\n", ctx->mc.sel_res_num);

    ctx->old_rotamer                       = ctx->cur_rotamers[ctx->mc.sel_res_num];
    ctx->cur_rotamers[ctx->mc.sel_res_num] = ctx->mc.sel_rotamer;

    for (i = 0; i < ctx->native_residue[ctx->mc.sel_res_num].ntorsions; i++) {
        a = integrator->sidechain_torsion[ctx->mc.sel_res_num][i][0];
        b = integrator->sidechain_torsion[ctx->mc.sel_res_num][i][1];
        c = integrator->sidechain_torsion[ctx->mc.sel_res_num][i][2];
        d = integrator->sidechain_torsion[ctx->mc.sel_res_num][i][3];

        if (integrator->USE_ROTAMERS)
            ctx->mc.delta_angle[i] =
                integrator->rotamer_angles[ctx->mc.sel_res_num].chis[sel_no_chi][i] +
                GaussianNum() *
                    sys->deviation_ang[ctx->native_residue[ctx->mc.sel_res_num].amino_num]
                                      [sel_no_chi][i] -
                ctx->native_residue[ctx->mc.sel_res_num].chi[i];
        else
            ctx->mc.delta_angle[i] = GaussianNum() * integrator->SIDECHAIN_NOISE;
        ctx->native_residue[ctx->mc.sel_res_num].tmpchi[i] =
            ctx->native_residue[ctx->mc.sel_res_num].chi[i] + ctx->mc.delta_angle[i];

        DoRotation(ctx, a, b, c, d, ctx->mc.delta_angle[i],
                   integrator->rotate_sidechain_natoms[ctx->mc.sel_res_num][i],
                   integrator->rotate_sidechain_atom[ctx->mc.sel_res_num][i]);
        for (j = 0; j < integrator->rotate_sidechain_natoms[ctx->mc.sel_res_num][i]; j++) {
            is_rotated[integrator->rotate_sidechain_atom[ctx->mc.sel_res_num][i][j]] = i + 1;
        }
    }

    ctx->all_rotated_natoms = integrator->rotate_sidechain_natoms[ctx->mc.sel_res_num][0];
    ctx->all_rotated_atoms  = integrator->rotate_sidechain_atom[ctx->mc.sel_res_num][0];

    return;
}

void SidechainMove(struct Context *ctx, struct Topology *top, struct MCIntegrator *integrator,
                   struct System *sys) {
    ctx->all_rotated_natoms = 0;
    ctx->total_pairs = ctx->total_pairs2 = 0;
    ctx->total_hbond_pairs               = 0;

    do {
        ctx->mc.sel_res_num = (int)(threefryrand() * top->nresidues);
    } while (ctx->native_residue[ctx->mc.sel_res_num].nrotamers <= 1 ||
             ctx->native_residue[ctx->mc.sel_res_num].amino_num == 14);

    MakeSidechainMove(ctx, top, integrator, sys);
    ctx->sidemovedone = 1;
    UpdateLattice(ctx, integrator, integrator->rotate_sidechain_natoms[ctx->mc.sel_res_num][0],
                  integrator->rotate_sidechain_atom[ctx->mc.sel_res_num][0]);
    NewDeltaContacts(ctx, sys, integrator->rotate_sidechain_natoms[ctx->mc.sel_res_num][0],
                     integrator->rotate_sidechain_atom[ctx->mc.sel_res_num][0],
                     integrator->sidechain_not_rotated[ctx->mc.sel_res_num][0]);

    return;
}
