#include "fold.h"
#include "define.h"
#include "constraint.h"
#include "contacts.h"
#include "energy.h"
#include "getrms.h"
#include "init.h"
#include "lattice_util.h"
#include "move.h"
#include "pdb_util.h"
#include "protein_util.h"
#include "hbonds.h"

/*
 * the_matrix is a 3D array of struct cell
 */
void SetupMatrixStuff(
    struct MCIntegrator *integrator,
    struct Topology *top,
    struct Context *ctx
) {
    int          i, j, k;
    struct cell *the_cell;

    for (i = 0; i < integrator->MATRIX_SIZE; i++)
        for (j = 0; j < integrator->MATRIX_SIZE; j++)
            for (k = 0; k < integrator->MATRIX_SIZE; k++)
                ctx->the_matrix[i][j][k].natoms = 0;

    for (i = 0; i < top->natoms; i++) {
        FindLatticeCoordinates(ctx, integrator, &ctx->native[i]);

        the_cell = &ctx->the_matrix[ctx->native[i].X][ctx->native[i].Y][ctx->native[i].Z];
        the_cell->atom_list[(the_cell->natoms)++] = i;
    }
}

void Fold(
    struct Context *ctx,
    struct Topology *top,
    struct Simulation *sim,
    struct System *sys,
    struct MCIntegrator *integrator
) {
    int i;
    int kkk;

    int   irep;
    int   itmp;
    int   ntmp;
    float etmp;
    int   partner;

    char  rmsd_filename[500];
    char  temp_filename[500];
    float E_pot_now        = 0.0;
    float E_tor_now        = 0.0;
    float E_sct_now        = 0.0;
    float E_aro_now        = 0.0;
    float E_hbond_now      = 0.0;
    float E_constraint_now = 0.0;  // AB
    ctx->disulfide_pairs        = (int *)calloc(sys->N_constraints, sizeof(int));

    clock_t t1, t1_f, t2, t2_f;
    /* rmsd by jyang */

    /* the following two moved to backbone.h */
    /* struct backbone struct_f1[MAXSEQUENCE], struct_f2[MAXSEQUENCE]; */

    /*
     * struct alignment contains int NFRAG as well as two arrays of
     * struct fragments of size MAXFRAG (=30).
     * struct fragments contains two ints: x1, x2
     */
    t1 = clock();
    struct alignment align;

    align.NFRAG           = 1;
    align.seqptr[1].x1    = 1;
    align.structptr[1].x1 = 1;
    align.seqptr[1].x2    = top->nresidues;
    align.structptr[1].x2 = top->nresidues;

    t1_f              = clock();
    double time_taken = ((double)t1_f - (double)t1) / CLOCKS_PER_SEC;
    printf("Initialization : %f seconds", time_taken);
    t2            = clock();
    ctx->mc_flags.init = !sim->NO_NEW_CLASHES;

    /* Save protein data into orig_native */
    for (i = 0; i < top->natoms; i++)
        CopyAtom(ctx->native[i], ctx->orig_native[i]);

    CenterProtein(&ctx->native, top->natoms);
    SetupMatrixStuff(integrator, top, ctx);

    Contacts(ctx, top, sys);
    fprintf(sim->STATUS, "INITIAL CLASHES\t%d\n", ctx->nclashes);
    if (!ctx->mc_flags.init) {
        fprintf(sim->STATUS, "\nturning off clashes:\n");
        TurnOffNativeClashes(ctx, top, sim, sys, 1);
    }
    fprintf(sim->STATUS, "\n");

    /* wmj */
    { /* Initialize contactstring */
        ctx->orig_contactstring = (short*)malloc(sizeof *ctx->orig_contactstring * top->nresidues * top->nresidues);
        fill_calpha_contact_string(
            top,
            ctx->native_residue, ctx->struct_native,
            ctx->orig_contactstring);
    }

    for (i = 0; i < top->nalign; ++i) {
        ctx->struct_f1[i + 1].CA = ctx->struct_native[ctx->struct_residue[top->map_to_struct[i]].CA].xyz;
        ctx->struct_f2[i + 1].CA = ctx->native[ctx->native_residue[top->map_to_seq[i]].CA].xyz;
    }
    ctx->native_rms = getrms(ctx->struct_f1, ctx->struct_f2, constraint_align);

    ResetEnergies(ctx, sim, sys, top, 0);  // Here we actually fill in the values for prev_E (which are, ironically, the current
             // energies, but these will indeed become previous once a move is made)
    ctx->rms_RMSDmin     = 100.;
    ctx->mcstep_RMSDmin  = 0;
    ctx->mcstep_Emin     = 0;
    ctx->Emin            = ctx->prev_E;
    ctx->Emin_pot        = ctx->prev_E_pot;
    ctx->Emin_hbond      = ctx->prev_E_hbond;
    ctx->Emin_tor        = ctx->prev_E_tor;
    ctx->Emin_sct        = ctx->prev_E_sct;
    ctx->Emin_aro        = ctx->prev_E_aro;
    ctx->Emin_constraint = ctx->prev_E_constraint;

    for (kkk = 0; kkk < sys->N_constraints; kkk++) {
        ctx->disulfide_pairs[kkk] = ctx->disulfide_pairs_attempt[kkk];
    }

    fprintf(sim->STATUS,
            "ENERGY = %.2f(clashes) + %.2f(rmsd) + %.2f(potential) + %.2f(Aro) + %.2f(hbond) + "
            "%.2f(tor) + %.2f(sct)\n\n",
            sys->weight_clash, sys->weight_rms, sys->weight_potential, ARO_WEIGHT, sys->weight_hbond, TOR_WEIGHT,
            SCT_WEIGHT);

    t2_f       = clock();
    time_taken = ((double)t2_f - (double)t2) / CLOCKS_PER_SEC;
    printf("Energy calculation : %f seconds", time_taken);

    for (i = 0; i < top->natoms; i++) {
        ctx->prev_native[i] = ctx->native[i];
        ctx->native_Emin[i] = ctx->prev_native[i];
    }
    for (i = 0; i < top->nresidues; ++i)
        ctx->cur_rotamers[i] = 0;

    sim->naccepted            = 0;
    sim->nrejected            = 0;
    sim->nothers              = 0;
    sim->n_sidechain_accepted = 0;

    ctx->natives = number_of_calpha_native_contacts(
        top,
        ctx->native_residue, ctx->native, ctx->orig_contactstring
    );
    fprintf(sim->STATUS,
            "         step #    energy contact rmsd   natives potnl    sctor      hbond       Aro  "
            " torsion accept reject others   temp setpoint E_constraint "
            "weighted_mean_constraint_distance\n---------------------------------------------------"
            "------------------------------------------\n");

    /* main folding loop */
    for (ctx->mcstep = 0; ctx->mcstep < sim->MC_STEPS; ctx->mcstep++) {
        /* AB: Turn off knowledge based moves at step MAX_CLUSTERSTEP */
        if (ctx->mcstep == integrator->MAX_CLUSTERSTEP) {
            integrator->USE_CLUSTER = 0;
            fprintf(sim->STATUS, "At MC step %10ld, USE_CLUSTER was set to %5.3f \n", ctx->mcstep,
                    integrator->USE_CLUSTER);
        } else if (ctx->mcstep == 0) {
            fprintf(sim->STATUS, "At MC step %10ld, USE_CLUSTER has value %5.3f \n", ctx->mcstep,
                    integrator->USE_CLUSTER);
        }

        /******** replica exchange ****************/
        if (ctx->mcstep % sim->MC_REPLICA_STEPS == 0) {
            // ierr=MPI_Barrier(mpi_world_comm);  //Added by AB to ensure synchronization--seems
            // like things were getting out of sync for some reason?

            for (i = 0; i < sim->nprocs; i++) {
                sim->replica_index[i]     = i;
                sim->attempted_replica[i] = 0;
            }

            sim->ierr = MPI_Allgather(&ctx->E, 1, MPI_FLOAT, sim->Enode, 1, MPI_FLOAT, sim->mpi_world_comm);
            sim->ierr = MPI_Allgather(&ctx->natives, 1, MPI_INT, sim->Nnode, 1, MPI_INT, sim->mpi_world_comm);

            if (sim->myrank == 0) {  // node 0 will mediate the exchange
                for (irep = 0; irep < sim->MAX_EXCHANGE; irep++) {
                    /* Choose replica pair for potential exchange: */
                    ctx->sel_num = (int)(threefryrand() * (sim->nprocs - 1)); /*Who initiates exchange?*/
                    /*AB: This was previously nprocs-2, but I saw no reason why second to last node
                    shouldn't be able to initiate exchange When we say (int)*random*(n_procs-1), we
                    keep in mind that (int) applies a floor function, so this will draw valeus
                    between 0 and nprocs-2 Moreover nprocs-2 corresponds to the second ot last node
                    due to zero indexing

                    Previously we were only drawing up to third to last node */

                    // Now, we need to choose the "partner" with whom sel_num exchanges
                    if (ctx->sel_num % sim->NODES_PER_TEMP ==
                        sim->NODES_PER_TEMP - 1) {  // the last setpoint for a given temperature
                        /*In this case, the exchange is initiated by the  lowest setpoint node in a
                        given temp It wouldn't't make sense for this node to exchange with the node
                        above it in index, since this would involve exchanging with  a partner who
                        has a very different setpoint than the initiator--the acceptance ratio woudl
                        be very low

                        Thus, we always exchange with the node that is one temperature step above it
                        (i.e. whose rank is NODES_PER_TEMP higher)
                        */
                        partner = ctx->sel_num + sim->NODES_PER_TEMP;
                    } else if (ctx->sel_num >= sim->nprocs - sim->NODES_PER_TEMP) {
                        /*This is saying that the exchange is initiated by somebody with the highest
                         * temperature, so obviously can't go up in temperature
                         */
                        partner = ctx->sel_num + 1;
                    } else {  // we either exchnage with the node that is one temperature above or
                              // one contact setpoint below, with 50% probability each
                        if (threefryrand() < 0.5) {
                            partner = ctx->sel_num + 1;
                        } else {
                            partner = ctx->sel_num + sim->NODES_PER_TEMP;
                        }
                    }

                    if (sim->attempted_replica[partner] == 1) {
                        /* Already some replica tried to exchange with partner */
                        continue;
                    }
                    sim->attempted_replica[partner] = 1;

                    // fprintf(STATUS, "irep : %d, sel_num : %d\n", irep, sel_num);
                    // fflush(STATUS);

                    /* Decide whether exchange will take place: */
                    ctx->delta_E = sim->Enode[partner] - sim->Enode[ctx->sel_num];
                    ctx->delta_T = 1.0 / sim->Tnode[partner] - 1.0 / sim->Tnode[ctx->sel_num];
                    ctx->delta_N =
                        sys->k_bias *
                            ((sim->Nnode[partner] - sim->Cnode[ctx->sel_num]) * (sim->Nnode[partner] - sim->Cnode[ctx->sel_num]) -
                             (sim->Nnode[ctx->sel_num] - sim->Cnode[ctx->sel_num]) *
                                 (sim->Nnode[ctx->sel_num] - sim->Cnode[ctx->sel_num])) +
                        sys->k_bias *
                            ((sim->Nnode[ctx->sel_num] - sim->Cnode[partner]) * (sim->Nnode[ctx->sel_num] - sim->Cnode[partner]) -
                             (sim->Nnode[partner] - sim->Cnode[partner]) * (sim->Nnode[partner] - sim->Cnode[partner]));
                    ctx->delta_all = ctx->delta_E * ctx->delta_T - ctx->delta_N;

                    if (ctx->delta_all >= 0 || threefryrand() < expf(ctx->delta_all)) {
                        itmp                   = sim->replica_index[ctx->sel_num];
                        sim->replica_index[ctx->sel_num] = sim->replica_index[partner];
                        sim->replica_index[partner] =
                            itmp;  // for instance, say 5 initiates exchange with 6...then
                                   // replica_index[5] now gets a value of 6, and vice versa
                        etmp           = sim->Enode[ctx->sel_num];
                        sim->Enode[ctx->sel_num] = sim->Enode[partner];
                        sim->Enode[partner] = etmp;
                        ntmp           = sim->Nnode[ctx->sel_num];
                        sim->Nnode[ctx->sel_num] = sim->Nnode[partner];
                        sim->Nnode[partner] = ntmp;
                        // accepted_replica[ctx->sel_num][ctx->partner]++;
                        /* Update count accordingly */
                        /* These counts are for exchanges upward (in temp or
                         * toward fewer native contacts). i.e. the counts of
                         * accept/reject are for the bottom partner. So you'll
                         * see 0 and 0 for the highest node (highest temperature
                         * and lowest native contact setpoint) */
                        sim->accepted_replica[ctx->sel_num]++;
                        // fprintf(STATUS, "Node %d successfully exchanged with Node %d, with a
                        // delta_all of %8.3f \n", sel_num, partner, delta_all);
                    } else {
                        // rejected_replica[sel_num][partner]++;
                        /* Update count accordingly */
                        sim->rejected_replica[ctx->sel_num]++;
                        // fprintf(STATUS, "Node %d UNSUCCESSFULLY exchanged with Node %d, with a
                        // delta_all of %8.3f \n", sel_num, partner, delta_all);
                    }
                }
            }

            /* At this point, all exchanges have been
             * decided. `replica_index` holds the new order of replicas,
             * e.g. `replica_index[i]` indicates which current replica will
             * go to replica i. */

            /* Let all nodes know; bcast from rank 0 */
            sim->ierr = MPI_Bcast(sim->replica_index, sim->nprocs, MPI_INT, 0, sim->mpi_world_comm);
            sim->ierr = MPI_Bcast(sim->Enode, sim->nprocs, MPI_FLOAT, 0, sim->mpi_world_comm);
            sim->ierr = MPI_Bcast(sim->Nnode, sim->nprocs, MPI_INT, 0, sim->mpi_world_comm);
            sim->ierr = MPI_Bcast(sim->accepted_replica, sim->nprocs, MPI_INT, 0, sim->mpi_world_comm);
            sim->ierr = MPI_Bcast(sim->rejected_replica, sim->nprocs, MPI_INT, 0, sim->mpi_world_comm);

            // fprintf(STATUS, "Replica Index\n");
            // for (i=0; i<nprocs; i++) fprintf(STATUS, "%5d %5d %8.3f\n", i, replica_index[i],
            // Enode[i]);

            /* Print output if it's time */
            if (ctx->mcstep % sim->MC_PRINT_STEPS == 0) {
                fprintf(sim->STATUS,
                        "RPLC %10ld E : %8.3f Natives : %d FROM %2d(T=%5.3f,setpoint=%d ) E : "
                        "%8.3f, Natives  :  %d, accepted : %5d, rejected : %5d  \n",
                        ctx->mcstep, sim->Enode[sim->myrank], sim->Nnode[sim->myrank], sim->replica_index[sim->myrank],
                        sim->Tnode[sim->replica_index[sim->myrank]], sim->Cnode[sim->replica_index[sim->myrank]],
                        sim->Enode[sim->myrank], sim->Nnode[sim->myrank], sim->accepted_replica[sim->myrank], sim->rejected_replica[sim->myrank]);
                // fprintf(STATUS, "RPLC %10ld E : %8.3f Natives : %d FROM %2d(T=%5.3f,setpoint=%d )
                // E : %8.3f, Natives  :  %d, accepted : %5d, rejected : %5d  \n", mcstep, E,
                // natives, replica_index[myrank], Tnode[replica_index[myrank]],
                // Cnode[replica_index[myrank]],Enode[myrank], Nnode[myrank],
                // accepted_replica[myrank][replica_index[myrank]],
                // rejected_replica[myrank][replica_index[myrank]]);
            }
            fflush(sim->STATUS);

            /* Actually do the exchange now */
            for (i = 0; i < top->natoms;
                 i++) {  // Temporary arrays to store all atom coordinates as they stood before any
                         // exchanges happened...these will be transferred later
                sim->buf_out.col(i) = ctx->native[i].xyz;
            }

            /* Keep in mind the following process is being carried out in
             * parallel. */

            /* Temporary holder of the coordinate replica number. i.e.,
             * the original replica number (from t = 0) of the coordinates
             * currently held by this replica is `current_replica`. */
            int receive_replicano = sim->current_replica;

            for (i = 0; i < sim->nprocs;
                 i++) {  // normally, replica_index[i] should equal i, unless an exchange occurred
                if (sim->replica_index[i] !=
                    i) {  // there is a discrepancy, indicating an exchange occurred...for instance,
                          // replica_index[5]=6 and replica_index[6]=5
                    if (sim->myrank == i) {  // For instance, if my rank is 6 and replica_index[6]=5,
                                        // then I received an exchange from node 5
                        sim->ierr = MPI_Recv(
                            sim->buf_in.data(), 3 * top->natoms, MPI_DOUBLE, sim->replica_index[i], (i + 2),
                            sim->mpi_world_comm,
                            &sim->mpi_status);  // now, receive the atomic coordinates from my partner
                        /* Also receive the "replica number" of the partner */
                        sim->ierr = MPI_Recv(&receive_replicano, 1, MPI_INT, sim->replica_index[i], (i + 3),
                                        sim->mpi_world_comm, &sim->mpi_status);
                    } else if (sim->myrank == sim->replica_index[i]) {  // For instance, if my rank is 5 and
                                                              // replica_index[6]=5, then this means
                                                              // I gave an exchange to node 6
                        sim->ierr = MPI_Send(sim->buf_out.data(), 3 * top->natoms, MPI_DOUBLE, i, (i + 2),
                                        sim->mpi_world_comm);  // give my coordinates to my partner
                        /* Send the "replica number" of this replica */
                        sim->ierr = MPI_Send(&sim->current_replica, 1, MPI_INT, i, (i + 3), sim->mpi_world_comm);
                    }
                }
            }
            sim->current_replica = receive_replicano;
            if (ctx->mcstep % sim->MC_PRINT_STEPS == 0) {
                fprintf(sim->STATUS, "REPNO %10ld %d \n", ctx->mcstep, sim->current_replica);
            }

            if (sim->replica_index[sim->myrank] !=
                sim->myrank) {  // if exchange occurred, transfer temporary data to real-deal atomic data
                           // structure, and recompute energies with the updated atomic info
                for (i = 0; i < top->natoms; i++) {
                    ctx->native[i].xyz = sim->buf_in.col(i);
                }
            }

            /* Re-center */
            //      if (!USE_ROTAMERS) {
            CenterProtein(&ctx->native, top->natoms);
            SetupMatrixStuff(integrator, top, ctx);
            for (i = 0; i < top->natoms; i++)
                ctx->prev_native[i] = ctx->native[i];
            //      }
            /* Re-set values that are susceptible to round-off error */
            Contacts(ctx, top, sys);
            if (!ctx->mc_flags.init)
                TurnOffNativeClashes(ctx, top, sim, sys, 0);
            ResetEnergies(ctx, sim, sys, top, 0);  // Recompute energies
            GetChi(top, ctx, integrator);

            ctx->natives = number_of_calpha_native_contacts(
                top, ctx->native_residue, ctx->native, ctx->orig_contactstring
            );
        } /************ END: replica exchange ********************/

        /* if backbone move was made and accepted we update native_rms */
        if (ctx->backbone_accepted == 1) {
            for (i = 0; i < top->nalign; ++i) {
                /* commented out because struct_native doesn't get altered */
                ctx->struct_f1[i + 1].CA = ctx->struct_native[ctx->struct_residue[top->map_to_struct[i]].CA].xyz;
                ctx->struct_f2[i + 1].CA = ctx->native[ctx->native_residue[top->map_to_seq[i]].CA].xyz;
            }
            ctx->native_rms = getrms(ctx->struct_f1, ctx->struct_f2, constraint_align);
        }

        /* Print Output */
        if (ctx->mcstep % sim->MC_PRINT_STEPS == 0) {
            // align_drms(native, native_residue, struct_native, struct_residue, map_to_seq,
            // map_to_struct, nalign, &native_rms);
            /* rmsd by jyang */
            for (i = 0; i < top->nalign; ++i) {
                /* commented out because struct_native doesn't get altered */
                ctx->struct_f1[i + 1].CA = ctx->struct_native[ctx->struct_residue[top->map_to_struct[i]].CA].xyz;
                ctx->struct_f2[i + 1].CA = ctx->native[ctx->native_residue[top->map_to_seq[i]].CA].xyz;
            }
            ctx->native_rms = getrms(ctx->struct_f1, ctx->struct_f2, constraint_align);
            TypeContacts(ctx, sys, top);
            E_pot_now   = FullAtomEnergy(sys, ctx);
            E_sct_now   = sctenergy(ctx, sys, top);
            E_tor_now   = torsionenergy(ctx, top, sys);
            E_aro_now   = aromaticenergy(ctx, sys, top);
            E_hbond_now = HydrogenBonds(top, ctx, sys);

            if (strcmp(sim->constraint_file, "None") != 0) {                                // AB
                E_constraint_now = Compute_constraint_energy(ctx, sys, ctx->native_residue, ctx->native);  // AB
            } else if (strcmp(sim->rmsd_constraint_file, "None") != 0) {
                E_constraint_now = Compute_frag_rmsd_constraint_energy(sys, ctx->struct_f1, ctx->struct_f2,
                                                                       constraint_align);  // AB
            } else {                                                                       // AB
                E_constraint_now = 0;                                                      // AB
            }  // AB

            fprintf(sim->STATUS,
                    "STEP %10ld  %8.2f %6d %5.2f %d %9.2f %9.2f %9.2f    %6.2f %8.2f %6.2f %6.2f "
                    "%6.2f %6.3f, %d %6.3f  \n",
                    ctx->mcstep, ctx->E, ctx->ncontacts, ctx->native_rms, ctx->natives, E_pot_now, E_sct_now, E_hbond_now,
                    E_aro_now, E_tor_now, 100 * (Float)sim->naccepted / (Float)sim->MC_PRINT_STEPS,
                    100 * (Float)sim->nrejected / (Float)sim->MC_PRINT_STEPS,
                    100 * (Float)sim->nothers / (Float)sim->MC_PRINT_STEPS, integrator->MC_TEMP,
                    sim->number_of_contacts_setpoint, E_constraint_now);
            fflush(sim->STATUS);

            sim->n_sidechain_accepted = 0;
            sim->naccepted            = 0;
            sim->nrejected            = 0;
            sim->nothers              = 0;
        }

        if (ctx->mcstep % 10000 == 0) {
            /* Re-center */
            CenterProtein(&ctx->native, top->natoms);
            SetupMatrixStuff(integrator, top, ctx);
            for (i = 0; i < top->natoms; i++)
                ctx->prev_native[i] = ctx->native[i];

            /* Re-set values that are susceptible to round-off error */
            Contacts(ctx, top, sys);
            if (!ctx->mc_flags.init)
                TurnOffNativeClashes(ctx, top, sim, sys, 0);
            ResetEnergies(ctx, sim, sys, top, ctx->mcstep);
            GetChi(top, ctx, integrator);
        }

        /* Output Structure */
        if (sim->PRINT_PDB) {
            if (ctx->mcstep % sim->MC_PDB_PRINT_STEPS == 0) {
                sprintf(temp_filename, "%s.%ld", sim->pdb_out_file, ctx->mcstep);
                //    if ((MC_TEMP < 0.18) || (mcstep > MC_STEPS - 100000))
                PrintPDB(sim, ctx, top, temp_filename);
            }
        }

        if (ctx->E < ctx->Emin) {
            ctx->mcstep_Emin     = ctx->mcstep;
            ctx->Emin            = ctx->E;
            ctx->Emin_pot        = ctx->E_pot;
            ctx->Emin_hbond      = ctx->E_hbond;
            ctx->Emin_tor        = ctx->E_tor;
            ctx->Emin_aro        = ctx->E_aro;
            ctx->Emin_constraint = ctx->E_constraint;  // AB
            ctx->rms_Emin        = ctx->native_rms;
            for (i = 0; i < top->natoms; i++)
                ctx->native_Emin[i] = ctx->native[i];
        }

        if (ctx->native_rms < ctx->rms_RMSDmin) {
            ctx->mcstep_RMSDmin       = ctx->mcstep;
            ctx->rms_RMSDmin          = ctx->native_rms;
            ctx->E_RMSDmin            = ctx->E;
            ctx->E_RMSDmin_pot        = ctx->E_pot;
            ctx->E_RMSDmin_hbond      = ctx->E_hbond;
            ctx->E_RMSDmin_tor        = ctx->E_tor;
            ctx->E_RMSDmin_sct        = ctx->E_sct;
            ctx->E_RMSDmin_aro        = ctx->E_aro;
            ctx->E_RMSDmin_constraint = ctx->E_constraint;  // AB
            for (i = 0; i < top->natoms; i++)
                ctx->native_RMSDmin[i] = ctx->native[i];
        }
        /* Make a move */
        MakeMove(
            ctx, sys, top, sim, integrator,
            integrator->STEP_SIZE, integrator->USE_GLOBAL_BB_MOVES
        );
    } /* end main folding loop */

    sprintf(temp_filename, "%s_Emin.pdb", sim->pdb_out_file);
    PrintPDB_Emin(sim, ctx, top, temp_filename);
    fprintf(sim->STATUS, "\nMC step at Emin: %10ld\n", ctx->mcstep_Emin);
    fprintf(sim->STATUS,
            "Emin:%8.2f  Emin_pot:%8.2f  E_hb:%7.2f  E_tor:%7.2f  E_sct:%7.2f E_aro:%7.2f "
            "E_constraint:%7.2f \n",
            ctx->Emin, ctx->Emin_pot, ctx->Emin_hbond, ctx->Emin_tor, ctx->Emin_sct, ctx->Emin_aro, ctx->Emin_constraint);
    fprintf(sim->STATUS, "rmsd at Emin:%8.2f  ", ctx->rms_Emin);
    fprintf(sim->STATUS, "Pdb file at Emin: %s\n", temp_filename);

    sprintf(rmsd_filename, "%s_RMSDmin.pdb", sim->pdb_out_file);
    PrintPDB_RMSDmin(sim, ctx, top, rmsd_filename);
    fprintf(sim->STATUS, "\nMC step at RMSDmin: %10ld\n", ctx->mcstep_RMSDmin);
    fprintf(sim->STATUS,
            "rms_Rmin:%8.2f  E_RMSDmin:%8.2f E_Rmin_pot:%8.2f E_Rhb:%8.2f E_Rtor:%8.2f "
            "E_Rsct:%8.2f E_Raro:%8.2f E_Rconstraint:%8.2f \n",
            ctx->rms_RMSDmin, ctx->E_RMSDmin, ctx->E_RMSDmin_pot, ctx->E_RMSDmin_hbond, ctx->E_RMSDmin_tor, ctx->E_RMSDmin_sct,
            ctx->E_RMSDmin_aro, ctx->E_RMSDmin_constraint);
    fprintf(sim->STATUS, "Pdb file at RMSDmin: %s\n", rmsd_filename);

    return;
}
