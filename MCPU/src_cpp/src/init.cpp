#define JSON_DIAGNOSTICS 1

#include "init.h"

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <cstring>
#include <format>

#include "constraint.h"
#include "contacts.h"
#include "define.h"  // Fixes float and constants
#include "energy.h"  // Fixes FullAtomEnergy, initialize_torsion, etc.
#include "fold.h"
#include "getrms.h"
#include "globals.h"       // Fixes nprocs, myrank, buf_in, amino_acids, etc.
#include "lattice_util.h"  // Fixes CenterProtein, get_template
#include "loop.h"
#include "pdb_util.h"
#include "protein_util.h"
#include "rng.h"  // Fixes threefryrand
#include "vector.h"
#include "align.h"
#include "hbonds.h"
#include "db/amino_acid_db.h"


#include <nlohmann/json.hpp>
#include <fstream>
#include <string>

using json = nlohmann::json;

short three[5] = {1, 3, 9, 27, 81};

/*================================================*/
/*           main initialization routine          */
/*================================================*/

void ReadTypesFile(struct Simulation *sim, struct Topology *top, struct System *sys) {
    int   type_num, i;
    char  res_name[4], atom_name[4];

    std::ifstream type_file(sim->atom_type_file);
    if (!type_file.is_open()) {
        fprintf(sim->STATUS, "ERROR: Can't open the file: %s!\n", sim->atom_type_file.c_str());
        exit(1);
    }
    i = -1;
    while (type_file.good()) {
        type_file >> atom_name >> res_name >> type_num;
        ++i;
        strcpy(top->atom_type_list[i].atom_name, atom_name);
        strcpy(top->atom_type_list[i].res_name, res_name);
        top->atom_type_list[i].type_num = type_num;
        if (!strcmp(res_name, "XXX")) {
            if (!strcmp(atom_name, "N"))
                top->bb_N_type = type_num;
            if (!strcmp(atom_name, "O"))
                top->bb_O_type = type_num;
            if (!strcmp(atom_name, "OXT"))
                top->bb_OXT_type = type_num;
        }
    }
    top->natom_type_list = i + 1;
    sys->MAX_TYPES       = top->atom_type_list[i].type_num + 1;
    type_file.close();
}

void SetupMuPotential(struct System *sys, struct Topology *top, struct Context *ctx,
                      struct Simulation *sim) {
    int   i, j, k, l;
    float n_nat, avg;
    int **nat_cons, **nnat_cons;

    nat_cons  = (int **)calloc(sys->MAX_TYPES, sizeof(int *));
    nnat_cons = (int **)calloc(sys->MAX_TYPES, sizeof(int *));
    for (i = 0; i < sys->MAX_TYPES; ++i) {
        nat_cons[i]  = (int *)calloc(sys->MAX_TYPES, sizeof(int));
        nnat_cons[i] = (int *)calloc(sys->MAX_TYPES, sizeof(int));
    }

    for (k = 0; k < top->natoms; ++k)
        for (l = k + 1; l < top->natoms; ++l) {
            if (ctx->data[k][l].contacts) {
                if (ctx->native[k].smogtype <= ctx->native[l].smogtype)
                    nat_cons[ctx->native[k].smogtype][ctx->native[l].smogtype]++;
                else
                    nat_cons[ctx->native[l].smogtype][ctx->native[k].smogtype]++;
            } else {
                if (ctx->native[k].smogtype <= ctx->native[l].smogtype)
                    nnat_cons[ctx->native[k].smogtype][ctx->native[l].smogtype]++;
                else
                    nnat_cons[ctx->native[l].smogtype][ctx->native[k].smogtype]++;
            }
        }

    avg   = 0;
    n_nat = 0;
    for (i = 0; i < sys->MAX_TYPES; ++i)
        for (j = i; j < sys->MAX_TYPES; ++j) {
            if (nat_cons[i][j] != 0) {
                avg += nnat_cons[i][j] / nat_cons[i][j];
                n_nat = n_nat + 1;
            }
        }
    sys->mu = (avg / (float)n_nat) / (1 + avg / (float)n_nat);
    fprintf(sim->STATUS, "mu: %f\n", sys->mu);

    for (i = 0; i < sys->MAX_TYPES; ++i)
        for (j = i; j < sys->MAX_TYPES; ++j)
            if ((nat_cons[i][j] != 0) || (nnat_cons[i][j] != 0)) {
                if (sys->mu == 1) {
                    if (nat_cons[i][j] != 0)
                        sys->potential[i][j] = -1;
                    else
                        sys->potential[i][j] = 1;
                } else if (sys->mu == 0) {
                    if (nnat_cons[i][j] != 0)
                        sys->potential[i][j] = 1;
                    else
                        sys->potential[i][j] = -1;
                } else /*AB: Here is the mu potential, ex. eq. (2) in Yang..Shakhnovich et. al PNAS
                          2007*/
                    sys->potential[i][j] =
                        ((1 - sys->mu) * nnat_cons[i][j] - sys->mu * nat_cons[i][j]) /
                        (sys->mu * nat_cons[i][j] + (1 - sys->mu) * nnat_cons[i][j]);
                sys->potential[j][i] = sys->potential[i][j];
            }

    for (i = 0; i < sys->MAX_TYPES; ++i) {
        free(nat_cons[i]);
        free(nnat_cons[i]);
    }
    free(nat_cons);
    free(nnat_cons);
}

void InitializeProtein(struct Context *ctx, struct Topology *top, struct Simulation *sim,
                       struct System *sys, struct MCIntegrator *integrator) {
    int i;
    /* Reset variables */
    /**See atom.h to learn about how atom structures are set up/
    fprintf(STATUS, "---MODEL---\n");
    fprintf(STATUS, "  file:\t\t%s   \n", native_file);  /*Set up arrays for native file, minimized
    file, and file with min RMSD*/
    ctx->native = (struct atom *)calloc(
        MAX_ATOMS, sizeof(struct atom)); /*Allocates memory for an array of atom_structures, each of
                                            which corresponds to data for one atom*/
    ctx->native_Emin    = (struct atom *)calloc(MAX_ATOMS, sizeof(struct atom));
    ctx->native_RMSDmin = (struct atom *)calloc(MAX_ATOMS, sizeof(struct atom));
    ctx->prev_native    = (struct atom *)calloc(MAX_ATOMS, sizeof(struct atom));
    ctx->orig_native    = (struct atom *)calloc(MAX_ATOMS, sizeof(struct atom));
    top->natoms         = 0;

    fprintf(sim->STATUS, "Initialized molecular data structures.\n");

    /* Initialize static data */

    // top->amino_acids = (struct amino *)calloc(20, sizeof(struct amino));
    top->amino_acids.resize(20);
    fprintf(sim->STATUS, "Made amino_acids structure\n");

    ReadTypesFile(sim, top, sys); /*Left off here*/
    fprintf(sim->STATUS, "Read Types File\n");

    /* Read data from native_file into native */
    ReadNative(sim, sys, integrator, top, sim->native_file.c_str(), ctx->native,
               &top->natoms);  // AB: This reads native_file into the structure native
    fprintf(sim->STATUS, "Read Native\n");

    /* Count nresidues */
    top->nresidues = 0;
    for (i = 0; i < top->natoms; i++)
        if (!strncmp(ctx->native[i].atomname, "CA", 2))
            top->nresidues++;

    if (top->nresidues != (ctx->native[top->natoms - 1].res_num + 1))
        fprintf(sim->STATUS, "FILE: %s   MISSING RESIDUES!!!\n", sim->native_file.c_str());

    fprintf(sim->STATUS, "  pdb length:\t\t%d\n  # of CA's:\t\t%d\n\n", top->nresidues,
            ctx->native[top->natoms - 1].res_num + 1);

    //sim->buf_in = Eigen::Matrix3Xd::Zero(3, top->natoms);
    //sim->buf_out = Eigen::Matrix3Xd::Zero(3, top->natoms);

    ReadAlignment(top, sys);
    fprintf(sim->STATUS, "Read Alignment\n");
    get_template(top, sim);
    fprintf(sim->STATUS, "Read Template Information\n");

    if (sys->USE_GO_POTENTIAL)
        sys->MAX_TYPES = top->natoms;

    SetRadii(sys);
    SetHardCore(sys, top, ctx);
    SetContactDistance(sys, top, ctx);

    /* Allocate potential-related structures */

    sys->potential = (float **)calloc(sys->MAX_TYPES, sizeof(float *));
    for (i = 0; i < sys->MAX_TYPES; i++) {
        sys->potential[i] = (float *)calloc(sys->MAX_TYPES, sizeof(float));
    }

    CenterProtein(&ctx->native, top->natoms);
    for (i = 0; i < top->natoms; i++)
        FindLatticeCoordinates(
            ctx, integrator,
            &ctx->native[i]); /* this also computes the integer value of the coordinates */

    /* Get residue info */

    ctx->native_residue = (struct residue *)calloc(
        top->nresidues,
        sizeof(struct residue));  // I believe native_residue stores the indices for the different
                                  // atoms (where they are located within native), and other info
    ctx->cur_rotamers = (int *)calloc(top->nresidues, sizeof(int));
    GetResidueInfo(ctx->native, ctx->native_residue, top->nresidues,
                   top->natoms);  // This adds the necessary info into native_residue

    /* Allocate memory for data structures */

    InitializeData(ctx, top, sim, sys);

    /* Determine phi-psi angles */

    GetPhiPsi(ctx->native, ctx->native_residue, top->nresidues);

    /* Set up correlation matrix */

    CheckCorrelation(ctx->data, ctx->native, ctx->native_residue, top->natoms, sim);

    /* Set up side-chain rotation structure */

    if (integrator->USE_SIDECHAINS) {
        /* Initialize sidechain rotation data */
        ReadSidechainTorsionData(sim, top);
        InitializeSidechainRotationData(top, integrator, ctx, sys, sim);
        /* initialize sidechain torsions */
        GetChi(top, ctx, integrator);
        if (integrator->USE_ROTAMERS) {
            ReadAvgChis(ctx, top, sim, integrator, sys);
            /* this will generate 1 side-chain move in the absence of backbone moves: */
            integrator->SIDECHAIN_MOVES = 1;
        }
    } else
        integrator->SIDECHAIN_MOVES = 0;

    /* initialize moved atoms data structures */

    ctx->ab = (struct pair *)calloc(top->natoms * top->natoms, sizeof(struct pair));
    ctx->cd = (struct pair *)calloc(top->natoms * top->natoms, sizeof(struct pair));

    initialize_torsion(top, sys, sim);
    initialize_sct(top, sys, sim);
    initialize_aromatic(top, sys, sim, ctx, integrator);
    initialize_secstr(sim, top);
    read_cluster(
        "/n/home01/kibumpark/pkg/dbfold2/MCPU/src_cpp/config_files/center30_state.csv",
        sys
    );
    std::cout << "DEBUG: Finished reading cluster data." << std::endl;
    if (sys->weight_hbond) {
        InitializeHydrogenBonding(top, ctx, sys, sim);
        std::cout << "DEBUG: Finished initializing hydrogen bonding data structures." << std::endl;
        ctx->hbond_pair = (struct pair *)calloc(top->natoms * top->natoms, sizeof(struct pair));
    }
    std::cout << "DEBUG: Finished initializing hydrogen bonding." << std::endl;

    /* Determine residue triplets */

    DetermineTriplets(top, integrator, sim, sys);

    std::cout << "DEBUG: Finished determining residue triplets." << std::endl;

    /* Initialize backbone rotation data */
    InitializeBackboneRotationData(top, integrator, ctx);

    std::cout << "DEBUG: Finished initializing backbone rotation data." << std::endl;

    /* Get contacts */

    Contacts(ctx, top, sys);

    std::cout << "DEBUG: Finished calculating contacts." << std::endl;

    /* Setup potential and related structures */

    SetupAlignmentStructure(ctx, top, sys, integrator, sim);

    std::cout << "DEBUG: Finished setting up alignment structure." << std::endl;

    if (sys->USE_GO_POTENTIAL) {
        sys->mu = 1;
        SetupAlignmentPotential(top, ctx, sys, sim);
    } else if (sim->READ_POTENTIAL)
        ReadPotential(sim, sys);
    else
        SetupMuPotential(sys, top, ctx, sim);

    TypeContacts(ctx, sys, top);
    ctx->native_E = FullAtomEnergy(sys, ctx);

    std::cout << "DEBUG: Finished setting up potential and calculating native energy." << std::endl;

    //  free(amino_acids);

    return;
}

void TurnOffNativeClashes(struct Context *ctx, struct Topology *top, struct Simulation *sim,
                          struct System *sys, int ReportBack) {
    int i, j;

    for (i = 0; i < top->natoms; i++)
        for (j = i + 1; j < top->natoms; j++)
            if (ctx->data[i][j].clashes) {
                if (ReportBack)
                    fprintf(sim->STATUS, "native clash\t%d %s - %s\t%d %s - %s \t %.3f %.3f\n", i,
                            ctx->native[i].res, ctx->native[i].atomname, j, ctx->native[j].res,
                            ctx->native[j].atomname,
                            sqrt(D2(ctx->native[i].xyz, ctx->native[j].xyz)),
                            sqrt(sys->hard_core[ctx->native[i].smogtype][ctx->native[j].smogtype]) /
                                100);
                ctx->data[i][j].clashes = ctx->data[j][i].clashes = 0;
                ctx->nclashes--;
            }
}

/*================================================*/
/*     routines for initializing static data      */
/*================================================*/

void SetHardCore(struct System *sys, struct Topology *top, struct Context *ctx) {
    int   i, j, k;
    float temp;
    float rad1, rad2;

    sys->hard_core = (long int **)calloc(sys->MAX_TYPES, sizeof(long int *));
    for (i = 0; i < sys->MAX_TYPES; i++)
        sys->hard_core[i] = (long int *)calloc(sys->MAX_TYPES, sizeof(long int));

    /* smogtype is always used for atom sizes */
    for (i = 0; i < sys->MAX_TYPES; i++)
        for (j = 0; j < sys->MAX_TYPES; j++) {
            if (!sys->USE_GO_POTENTIAL) {
                for (k = 0; k < top->natom_type_list; ++k)
                    if (i == top->atom_type_list[k].type_num)
                        break;
                rad1 = sys->radii[TypeAtom(top->atom_type_list[k].atom_name,
                                           top->atom_type_list[k].res_name)];
                for (k = 0; k < top->natom_type_list; ++k)
                    if (j == top->atom_type_list[k].type_num)
                        break;
                rad2 = sys->radii[TypeAtom(top->atom_type_list[k].atom_name,
                                           top->atom_type_list[k].res_name)];
            } else {
                rad1 = sys->radii[ctx->native[i].atomtype];
                rad2 = sys->radii[ctx->native[j].atomtype];
            }
            temp                 = sys->ALPHA * (rad1 + rad2);
            sys->hard_core[i][j] = (long int)(temp * temp * INT_PRECISION * INT_PRECISION);
        }

    return;
}

void SetContactDistance(struct System *sys, struct Topology *top, struct Context *ctx) {
    float rad1, rad2;
    float temp;
    float temp00;
    int   i, j, k;

    sys->contact_distance = (struct cutoff **)calloc(sys->MAX_TYPES, sizeof(struct cutoff *));
    for (i = 0; i < sys->MAX_TYPES; i++)
        sys->contact_distance[i] = (struct cutoff *)calloc(sys->MAX_TYPES, sizeof(struct cutoff));

    for (i = 0; i < sys->MAX_TYPES; i++)
        for (j = 0; j < sys->MAX_TYPES; j++) {
            if (!sys->USE_GO_POTENTIAL) {
                for (k = 0; k < top->natom_type_list; ++k)
                    if (i == top->atom_type_list[k].type_num)
                        break;
                rad1 = sys->radii[TypeAtom(top->atom_type_list[k].atom_name,
                                           top->atom_type_list[k].res_name)];
                for (k = 0; k < top->natom_type_list; ++k)
                    if (j == top->atom_type_list[k].type_num)
                        break;
                rad2 = sys->radii[TypeAtom(top->atom_type_list[k].atom_name,
                                           top->atom_type_list[k].res_name)];
            } else {
                rad1 = sys->radii[ctx->native[i].atomtype];
                rad2 = sys->radii[ctx->native[j].atomtype];
            }
            temp                          = sys->LAMBDA * sys->ALPHA * (rad1 + rad2);
            sys->contact_distance[i][j].b = (long int)(temp * temp * INT_PRECISION * INT_PRECISION);
            sys->contact_distance[i][j].a = 0;
            sys->beta                     = 0.;
            temp00                        = sys->beta * (rad1 + rad2);
            sys->contact_distance[i][j].a =
                (long int)(temp00 * temp00 * INT_PRECISION * INT_PRECISION);
        }

    /* sets up proper distances for h-bonding */
    if (!sys->USE_GO_POTENTIAL) {
        /* turns off all backbone-sidechain contacts */
        for (j = top->bb_N_type; j <= top->bb_OXT_type; ++j)
            for (i = 0; i < sys->MAX_TYPES; ++i)
                if ((i < top->bb_N_type) || (i > top->bb_OXT_type)) {
                    sys->contact_distance[i][j].a = sys->contact_distance[j][i].a = 0;
                    sys->contact_distance[i][j].b = sys->contact_distance[j][i].b = 0;
                }
    }
    return;
}

void SetRadii(struct System *sys) {
    sys->radii = (float *)calloc(13, sizeof(float));

    /* carbons */
    sys->radii[0] = 1.61;
    sys->radii[1] = 1.76;
    sys->radii[2] = sys->radii[3] = sys->radii[4] = 1.88;
    /* nitrogens */
    sys->radii[5] = sys->radii[6] = sys->radii[7] = sys->radii[8] = 1.64;
    /* oxygens */
    sys->radii[9]  = 1.42;
    sys->radii[10] = 1.46;
    /* sulfurs */
    sys->radii[11] = sys->radii[12] = 1.77;

    return;
}

/*=============================================*/
/*   routines for initializing and reading     */
/*               the native pdb                */
/*=============================================*/

void ReadNative(struct Simulation *sim,
    struct System *sys, struct MCIntegrator *integrator, struct Topology *top,
    const char *file_name, struct atom *protein, int *Natoms) {
    FILE *the_file;
    char  line[250];

    *Natoms = 0;
    if ((the_file = fopen(file_name, "r")) == NULL) {
        fprintf(sim->STATUS, "ERROR: Can't open the file: %s!\n", file_name);
        exit(1);
    }
    while (fgets(line, 100, the_file) != NULL)
        ParsePDBLine(sim, sys, integrator, top, line, protein, Natoms);
    fclose(the_file);

    return;
}

void GetResidueInfo(struct atom *Chain, struct residue *Residue, int Nres, int Natoms) {
    int i, j;

    for (i = 0; i < Natoms; i++)
        for (j = 0; j < 40; j++)
            Residue[Chain[i].res_num].atomnumber[j] = -1;

    for (i = 0; i < Natoms; i++) {
        j                                       = MatchAtomname(Chain[i].atomname);
        Residue[Chain[i].res_num].atomnumber[j] = i;
        if (!strncmp(Chain[i].atomname, "CA", 2)) {
            Residue[Chain[i].res_num].CA = i;
            strcpy(Residue[Chain[i].res_num].res, Chain[i].res);
            Residue[Chain[i].res_num].amino_num   = GetAminoNumber(Chain[i].res);
            Residue[Chain[i].res_num].psi         = 0;
            Residue[Chain[i].res_num].phi         = 0;
            Residue[Chain[i].res_num].chi[0]      = 0;
            Residue[Chain[i].res_num].chi[1]      = 0;
            Residue[Chain[i].res_num].chi[2]      = 0;
            Residue[Chain[i].res_num].chi[3]      = 0;
            Residue[Chain[i].res_num].is_core     = Chain[i].is_core;
            Residue[Chain[i].res_num].is_designed = Chain[i].is_designed;
        } else if (!strcmp(Chain[i].atomname, "N"))
            Residue[Chain[i].res_num].N = i;
        else if (!strcmp(Chain[i].atomname, "C"))
            Residue[Chain[i].res_num].C = i;
        else if (!strcmp(Chain[i].atomname, "O"))
            Residue[Chain[i].res_num].O = i;
        else if (!strcmp(Chain[i].atomname, "CB"))
            Residue[Chain[i].res_num].CB = i;
        else if (!strcmp(Chain[i].atomname, "CG"))
            Residue[Chain[i].res_num].CG = i;
        else if (!strcmp(Chain[i].atomname, "CE1"))
            Residue[Chain[i].res_num].CE1 = i;
        else if (!strcmp(Chain[i].atomname, "CE2"))
            Residue[Chain[i].res_num].CE2 = i;
        else if (!strcmp(Chain[i].atomname, "CZ2"))
            Residue[Chain[i].res_num].CZ2 = i;
        else if (!strcmp(Chain[i].atomname, "CZ3"))
            Residue[Chain[i].res_num].CZ3 = i;
    }

    /* Handle CB atoms for glycine */

    for (i = 0; i < Nres; i++) {
        if (!strcmp(Residue[i].res, "GLY")) {
            Residue[i].CB = -999;
        }
        /*  AddCB(native, native_residue[i]); */
    }

    return;
}

void GetPhiPsi(struct atom *Chain, struct residue *Residue,
               int Nres) { /* Phi/Psi angles are stored in radians */
    int i;

    for (i = 0; i < Nres; i++) {
        if (i != 0)
            Residue[i].phi = PI / 180.0 * Phi(Residue[i], Residue[i - 1], Chain);
        else
            Residue[i].phi = -999;
        if (i != Nres - 1)
            Residue[i].psi = PI / 180.0 * Psi(Residue[i], Residue[i + 1], Chain);
        else
            Residue[i].psi = -999;
    }

    return;
}

void GetChi(struct Topology *top, struct Context *ctx,
            struct MCIntegrator *integrator) { /* Chi angles are stored in radians */
    int i, j;
    /* this routine also resets the native chis to the current values */

    for (i = 0; i < top->nresidues; i++)
        for (j = 0; j < top->amino_acids[ctx->native_residue[i].amino_num].ntorsions; j++) {
            ctx->native_residue[i].chi[j] =
                PI / 180.0 *
                CalculateTorsion(ctx->native, integrator->sidechain_torsion[i][j][0],
                                 integrator->sidechain_torsion[i][j][1],
                                 integrator->sidechain_torsion[i][j][2],
                                 integrator->sidechain_torsion[i][j][3], 0);
            ctx->native_residue[i].native_chi[j] = ctx->native_residue[i].chi[j];
            ctx->native_residue[i].tmpchi[j]     = ctx->native_residue[i].chi[j];
        }

    return;
}

void ReadAvgChis(struct Context *ctx, struct Topology *top, struct Simulation *sim,
                 struct MCIntegrator *integrator, struct System *sys) {
    int   i;
    char  name[4], line[200];
    float X, Y, Z, W;
    float value;
    float sX, sY, sZ, sW;

    integrator->rotamer_angles = (struct angles*)calloc(top->nresidues, sizeof(struct angles));

    if ((sim->DATA = fopen(sim->rotamer_data_file.c_str(), "r")) == NULL) {
        fprintf(sim->STATUS, "ERROR: Can't open the file: %s!\n", sim->rotamer_data_file);
        exit(1);
    }
    while (fgets(line, 150, sim->DATA) != NULL) {
        sscanf(line, "%s %*d%*d%*d%*d %*d%*d %f%*f%*f%*f  %f%f %f%f %f%f %f%f", name, &value, &X,
               &sX, &Y, &sY, &Z, &sZ, &W, &sW);
        for (i = 0; i < top->nresidues; ++i) {
            if (strcmp(name, ctx->native_residue[i].res) == 0) {
                integrator->rotamer_angles[i].chis[top->no_chi_list[GetAminoNumber(name)]][0] = X;
                integrator->rotamer_angles[i].chis[top->no_chi_list[GetAminoNumber(name)]][1] = Y;
                integrator->rotamer_angles[i].chis[top->no_chi_list[GetAminoNumber(name)]][2] = Z;
                integrator->rotamer_angles[i].chis[top->no_chi_list[GetAminoNumber(name)]][3] = W;
            }
        }
        sys->deviation_ang[GetAminoNumber(name)][top->no_chi_list[GetAminoNumber(name)]][0] = sX;
        sys->deviation_ang[GetAminoNumber(name)][top->no_chi_list[GetAminoNumber(name)]][1] = sY;
        sys->deviation_ang[GetAminoNumber(name)][top->no_chi_list[GetAminoNumber(name)]][2] = sZ;
        sys->deviation_ang[GetAminoNumber(name)][top->no_chi_list[GetAminoNumber(name)]][3] = sW;
        sys->prob_ang[GetAminoNumber(name)][top->no_chi_list[GetAminoNumber(name)]]         = value;
        top->no_chi_list[GetAminoNumber(name)]++;
    }
    fclose(sim->DATA);
}

void ReadSidechainTorsionData(struct Simulation *sim, struct Topology *top) {
    // 1. Load the library from JSON
    auto aa_library = mcpu::db::LoadTopology("/n/home01/kibumpark/pkg/dbfold2/MCPU/src_cpp/config_files/amino_torsion.json");
    if (aa_library.empty()) {
        fprintf(sim->STATUS, "ERROR: Can't load amino acid torsion library!\n");
        exit(1);
    }

    // Ensure the topology vector is sized for the 20 standard amino acids
    if (top->amino_acids.size() < 20) {
        top->amino_acids.resize(20);
    }

    for (const auto& aa : aa_library) {
        // 2. Identify the index (using aa.name from our JSON)
        int aa_num = GetAminoNumber(const_cast<char*>(aa.name.c_str()));
        
        if (aa_num < 0 || aa_num >= 20) {
            fprintf(sim->STATUS, "ERROR: Invalid amino acid in torsion library: %s\n", aa.name.c_str());
            continue; // Or exit(1) depending on strictness
        }

        auto& target = top->amino_acids[aa_num];

        // 3. Copy Metadata (Safe copy to std::array<char, N>)
        std::fill(target.name.begin(), target.name.end(), '\0');
        std::strncpy(target.name.data(), aa.name.c_str(), 3);
        
        std::fill(target.symbol.begin(), target.symbol.end(), '\0');
        std::strncpy(target.symbol.data(), aa.symbol.c_str(), 1);

        target.ntorsions = aa.ntorsions;

        // 4. Rotamer Logic (MCPU Specific)
        if (aa.name == "PRO") {
            target.nrotamers = 2;
        } else if (aa.name == "TYR" || aa.name == "HIS" || aa.name == "PHE") {
            target.nrotamers = 6;
        } else if (aa.name == "GLY") {
            target.nrotamers = 0;
        } else {
            // 'three' is the precomputed powers of 3 array (3^ntorsions)
            target.nrotamers = (int)three[aa.ntorsions];
        }

        // 5. Map Torsion Definitions and Rotation Atoms
        for (int t_idx = 0; t_idx < aa.ntorsions && t_idx < 4; ++t_idx) {
            const auto& t_src = aa.torsions[t_idx];

            // Define the 4 atoms forming the torsion (N, CA, CB, CG, etc.)
            for (int k = 0; k < 4 && k < t_src.atoms.size(); ++k) {
                std::strncpy(target.torsion[t_idx][k].data(), t_src.atoms[k].c_str(), 3);
                target.torsion[t_idx][k][3] = '\0';
            }

            // Define the atoms that move when this torsion rotates
            target.rotate_natoms[t_idx] = (int)t_src.affected_atoms.size();
            for (int k = 0; k < target.rotate_natoms[t_idx] && k < 10; ++k) {
                std::strncpy(target.rotate_atom[t_idx][k].data(), t_src.affected_atoms[k].c_str(), 3);
                target.rotate_atom[t_idx][k][3] = '\0';
            }
        }
    }

    fprintf(sim->STATUS, "DEBUG: Sidechain torsion data synchronized for %zu residues.\n", aa_library.size());
}

/*=============================================================*/
/*                Read full atom potential                     */
/*=============================================================*/

void ReadPotential(struct Simulation *sim, struct System *sys) {
    FILE *pot_file;
    int   i, j;
    float val;
    /* read max types potential */

    if ((pot_file = fopen(sim->potential_file.c_str(), "r")) == NULL) {
        fprintf(sim->STATUS, "ERROR: Can't open the file: %s!\n", sim->potential_file);
        exit(1);
    }
    while (fscanf(pot_file, "%d %d %f", &i, &j, &val) != EOF) {
        sys->potential[i][j] = sys->potential[j][i] = val;
    }

    fclose(pot_file);
}

/*=============================================================*/
/*                   allocate data structures                  */
/*=============================================================*/

void InitializeData(struct Context *ctx, struct Topology *top, struct Simulation *sim,
                    struct System *sys) {
    int i;
#if DEBUG
    sim->debug_contacts = (unsigned char **)calloc(top->natoms, sizeof(unsigned char *));
    for (i = 0; i < top->natoms; i++)
        sim->debug_contacts[i] = (unsigned char *)calloc(top->natoms, sizeof(unsigned char));
    sim->debug_dcontacts = (unsigned char **)calloc(top->natoms, sizeof(unsigned char *));
    for (i = 0; i < top->natoms; i++)
        sim->debug_dcontacts[i] = (unsigned char *)calloc(top->natoms, sizeof(unsigned char));
    sim->debug_clashes = (unsigned char **)calloc(top->natoms, sizeof(unsigned char *));
    for (i = 0; i < top->natoms; i++)
        sim->debug_clashes[i] = (unsigned char *)calloc(top->natoms, sizeof(unsigned char));
#endif

    ctx->data = (struct contact_data **)calloc(top->natoms, sizeof(struct contact_data *));
    for (i = 0; i < top->natoms; i++)
        ctx->data[i] = (struct contact_data *)calloc(top->natoms, sizeof(struct contact_data));

    ctx->type_contacts = (short **)calloc(sys->MAX_TYPES, sizeof(short *));
    for (i = 0; i < sys->MAX_TYPES; i++)
        ctx->type_contacts[i] = (short *)calloc(sys->MAX_TYPES, sizeof(short));

    is_rotated = (unsigned char *)calloc(top->natoms, sizeof(unsigned char));
    for (i = 0; i < top->natoms; i++)
        is_rotated[i] = 0;

    return;
}

/*=============================================================*/
/*           routines for initializing the move set            */
/*=============================================================*/

void DetermineTriplets(struct Topology *top, struct MCIntegrator *integrator,
                       struct Simulation *sim, struct System *sys) {
    int i, j, k;

    integrator->residue_triplets =
        (struct triplet *)calloc((top->nresidues - 5) * 81, sizeof(struct triplet));
    /* this is way more memory than needed, given that the exact number of triplets is 16n - 55 */
    integrator->total_triplets = 0;
    for (i = 0; i < top->nresidues; i++) {
        integrator->residue_triplets[integrator->total_triplets].a   = i;
        integrator->residue_triplets[integrator->total_triplets].b   = -1;
        integrator->residue_triplets[integrator->total_triplets++].c = -1;
    }
    integrator->TOTAL_SINGLE_LOOP_MOVES = integrator->total_triplets;

    for (i = 0; i < top->nresidues; i++)
        for (j = i + 1; j < top->nresidues; j++)
            if (j < i + 6) {
                integrator->residue_triplets[integrator->total_triplets].a   = i;
                integrator->residue_triplets[integrator->total_triplets].b   = j;
                integrator->residue_triplets[integrator->total_triplets++].c = -1;
            }
    integrator->TOTAL_DOUBLE_LOOP_MOVES =
        integrator->total_triplets - integrator->TOTAL_SINGLE_LOOP_MOVES;

    for (i = 0; i < top->nresidues; i++)
        for (j = i + 1; j < top->nresidues; j++)
            for (k = j + 1; k < top->nresidues; k++)
                if (k < i + 6 && j < i + 6) {
                    integrator->residue_triplets[integrator->total_triplets].a   = i;
                    integrator->residue_triplets[integrator->total_triplets].b   = j;
                    integrator->residue_triplets[integrator->total_triplets++].c = k;
                }
    integrator->TOTAL_TRIPLE_LOOP_MOVES = integrator->total_triplets -
                                          integrator->TOTAL_SINGLE_LOOP_MOVES -
                                          integrator->TOTAL_DOUBLE_LOOP_MOVES;

    return;
}

void InitializeBackboneRotationData(struct Topology *top, struct MCIntegrator *integrator,
                                    struct Context *ctx) {
    int i, j, k;

    /* rotate_atom[0=psi, 1=phi][which residue][list of atoms] */

    integrator->yang_rotated_atoms = (short *)calloc(top->natoms, sizeof(short));
    integrator->yang_not_rotated   = (char *)calloc(top->natoms, sizeof(char));

    integrator->rotate_natoms = (short **)calloc(2, sizeof(short *));
    integrator->rotate_atom   = (short ***)calloc(2, sizeof(short **));
    integrator->not_rotated   = (char ***)calloc(2, sizeof(char **));
    for (i = 0; i < 2; i++) {
        integrator->rotate_atom[i]   = (short **)calloc(top->nresidues, sizeof(short *));
        integrator->not_rotated[i]   = (char **)calloc(top->nresidues, sizeof(char *));
        integrator->rotate_natoms[i] = (short *)calloc(top->nresidues, sizeof(short));
        for (j = 0; j < top->nresidues; j++) {
            integrator->rotate_atom[i][j] = (short *)calloc(top->natoms, sizeof(short));
            integrator->not_rotated[i][j] = (char *)calloc(top->natoms, sizeof(char));
        }
    }

    /* rotate the short end of the chain for each residue */
    /* and determine which atoms were rotated, for either phi or psi rotation */

    for (i = 0; i < top->nresidues; i++) {
        integrator->rotate_natoms[0][i] = 0;
        integrator->rotate_natoms[1][i] = 0;

        if (i > top->nresidues / 2.0) {
            for (j = 0; j < top->natoms; j++)
                if (ctx->native[j].res_num > i) {
                    integrator->rotate_atom[0][i][integrator->rotate_natoms[0][i]++] = j;
                    integrator->rotate_atom[1][i][integrator->rotate_natoms[1][i]++] = j;
                } else if (ctx->native[j].res_num == i) {
                    if (j != ctx->native_residue[i].N && j != ctx->native_residue[i].CA)
                        integrator->rotate_atom[1][i][integrator->rotate_natoms[1][i]++] = j;
                    if (j == ctx->native_residue[i].O)
                        integrator->rotate_atom[0][i][integrator->rotate_natoms[0][i]++] = j;
                }
        } else {
            for (j = 0; j < top->natoms; j++)
                if (ctx->native[j].res_num < i) {
                    integrator->rotate_atom[0][i][integrator->rotate_natoms[0][i]++] = j;
                    integrator->rotate_atom[1][i][integrator->rotate_natoms[1][i]++] = j;
                } else if (ctx->native[j].res_num == i) {
                    if (j != ctx->native_residue[i].O && j != ctx->native_residue[i].C &&
                        j != ctx->native_residue[i].CA)
                        integrator->rotate_atom[0][i][integrator->rotate_natoms[0][i]++] = j;
                }
        }
    }

    /* set up the not_rotated array, with 1's at every unrotated atom */

    for (i = 0; i < top->nresidues; i++) {
        for (j = 0; j < integrator->rotate_natoms[0][i]; j++)
            integrator->not_rotated[0][i][integrator->rotate_atom[0][i][j]] = 1;
        for (j = 0; j < integrator->rotate_natoms[1][i]; j++)
            integrator->not_rotated[1][i][integrator->rotate_atom[1][i][j]] = 1;
        for (j = 0; j < top->natoms; j++) {
            integrator->not_rotated[0][i][j] = !integrator->not_rotated[0][i][j];
            integrator->not_rotated[1][i][j] = !integrator->not_rotated[1][i][j];
        }
    }

    /* store atoms rotated by loop moves */

    integrator->loop_rotate_natoms = (short **)calloc(integrator->total_triplets, sizeof(short *));
    integrator->loop_int_rotate_natoms =
        (short **)calloc(integrator->total_triplets, sizeof(short *));
    integrator->loop_rotate_atoms = (short ***)calloc(integrator->total_triplets, sizeof(short **));
    integrator->loop_int_rotate_atoms =
        (short ***)calloc(integrator->total_triplets, sizeof(short **));
    integrator->loop_not_rotated = (char ***)calloc(integrator->total_triplets, sizeof(char **));
    for (i = 0; i < integrator->total_triplets; i++) {
        /* there are up to 6 bonds to rotate for each loop move */
        integrator->loop_rotate_natoms[i]     = (short *)calloc(6, sizeof(short));
        integrator->loop_int_rotate_natoms[i] = (short *)calloc(6, sizeof(short));
        integrator->loop_rotate_atoms[i]      = (short **)calloc(6, sizeof(short *));
        integrator->loop_int_rotate_atoms[i]  = (short **)calloc(6, sizeof(short *));
        integrator->loop_not_rotated[i]       = (char **)calloc(6, sizeof(char *));

        /* all loop moves have at least 2 bonds */
        integrator->loop_rotate_atoms[i][0]     = (short *)calloc(top->natoms, sizeof(short));
        integrator->loop_rotate_atoms[i][1]     = (short *)calloc(top->natoms, sizeof(short));
        integrator->loop_int_rotate_atoms[i][0] = (short *)calloc(top->natoms, sizeof(short));
        integrator->loop_int_rotate_atoms[i][1] = (short *)calloc(top->natoms, sizeof(short));
        integrator->loop_not_rotated[i][0]      = (char *)calloc(top->natoms, sizeof(char));
        integrator->loop_not_rotated[i][1]      = (char *)calloc(top->natoms, sizeof(char));

        if (integrator->residue_triplets[i].b >= 0) {
            /* double and triple loop moves */
            integrator->loop_rotate_atoms[i][2]     = (short *)calloc(top->natoms, sizeof(short));
            integrator->loop_rotate_atoms[i][3]     = (short *)calloc(top->natoms, sizeof(short));
            integrator->loop_int_rotate_atoms[i][2] = (short *)calloc(top->natoms, sizeof(short));
            integrator->loop_int_rotate_atoms[i][3] = (short *)calloc(top->natoms, sizeof(short));
            integrator->loop_not_rotated[i][2]      = (char *)calloc(top->natoms, sizeof(char));
            integrator->loop_not_rotated[i][3]      = (char *)calloc(top->natoms, sizeof(char));
        } else {
            /* single loop moves */
            integrator->loop_rotate_atoms[i][2] =
                (short *)calloc(1, sizeof(short)); /* why are these allocated? */
            integrator->loop_rotate_atoms[i][3]     = (short *)calloc(1, sizeof(short));
            integrator->loop_int_rotate_atoms[i][2] = (short *)calloc(1, sizeof(short));
            integrator->loop_int_rotate_atoms[i][3] = (short *)calloc(1, sizeof(short));
            integrator->loop_not_rotated[i][2]      = (char *)calloc(top->natoms, sizeof(char));
            integrator->loop_not_rotated[i][3]      = (char *)calloc(top->natoms, sizeof(char));
        }

        if (integrator->residue_triplets[i].c >= 0) {
            /* triple loop moves */
            integrator->loop_rotate_atoms[i][4]     = (short *)calloc(top->natoms, sizeof(short));
            integrator->loop_rotate_atoms[i][5]     = (short *)calloc(top->natoms, sizeof(short));
            integrator->loop_int_rotate_atoms[i][4] = (short *)calloc(top->natoms, sizeof(short));
            integrator->loop_int_rotate_atoms[i][5] = (short *)calloc(top->natoms, sizeof(short));
            integrator->loop_not_rotated[i][4]      = (char *)calloc(top->natoms, sizeof(char));
            integrator->loop_not_rotated[i][5]      = (char *)calloc(top->natoms, sizeof(char));
        } else {
            /* single and double loop moves */
            integrator->loop_rotate_atoms[i][4]     = (short *)calloc(1, sizeof(short));
            integrator->loop_rotate_atoms[i][5]     = (short *)calloc(1, sizeof(short));
            integrator->loop_int_rotate_atoms[i][4] = (short *)calloc(1, sizeof(short));
            integrator->loop_int_rotate_atoms[i][5] = (short *)calloc(1, sizeof(short));
            integrator->loop_not_rotated[i][4]      = (char *)calloc(top->natoms, sizeof(char));
            integrator->loop_not_rotated[i][5]      = (char *)calloc(top->natoms, sizeof(char));
        }
    }

    for (i = 0; i < integrator->total_triplets; i++) {
        if (integrator->residue_triplets[i].a > top->nresidues / 2.0) { /* phi then psi */

            integrator->loop_rotate_natoms[i][0] =
                integrator->rotate_natoms[1][integrator->residue_triplets[i].a]; /* phi */
            for (j = 0; j < integrator->loop_rotate_natoms[i][0]; j++)
                integrator->loop_rotate_atoms[i][0][j] =
                    integrator->rotate_atom[1][integrator->residue_triplets[i].a][j];
            integrator->loop_rotate_natoms[i][1] =
                integrator->rotate_natoms[0][integrator->residue_triplets[i].a]; /* psi */
            for (j = 0; j < integrator->loop_rotate_natoms[i][1]; j++)
                integrator->loop_rotate_atoms[i][1][j] =
                    integrator->rotate_atom[0][integrator->residue_triplets[i].a][j];

            if (integrator->residue_triplets[i].b >= 0) {
                integrator->loop_rotate_natoms[i][2] =
                    integrator->rotate_natoms[1][integrator->residue_triplets[i].b]; /* phi */
                for (j = 0; j < integrator->loop_rotate_natoms[i][2]; j++)
                    integrator->loop_rotate_atoms[i][2][j] =
                        integrator->rotate_atom[1][integrator->residue_triplets[i].b][j];
                integrator->loop_rotate_natoms[i][3] =
                    integrator->rotate_natoms[0][integrator->residue_triplets[i].b]; /* psi */
                for (j = 0; j < integrator->loop_rotate_natoms[i][3]; j++)
                    integrator->loop_rotate_atoms[i][3][j] =
                        integrator->rotate_atom[0][integrator->residue_triplets[i].b][j];
            } else {
                integrator->loop_rotate_natoms[i][2] = 0;
                integrator->loop_rotate_natoms[i][3] = 0;
            }
            if (integrator->residue_triplets[i].c >= 0) {
                integrator->loop_rotate_natoms[i][4] =
                    integrator->rotate_natoms[1][integrator->residue_triplets[i].c]; /* phi */
                for (j = 0; j < integrator->loop_rotate_natoms[i][4]; j++)
                    integrator->loop_rotate_atoms[i][4][j] =
                        integrator->rotate_atom[1][integrator->residue_triplets[i].c][j];
                integrator->loop_rotate_natoms[i][5] =
                    integrator->rotate_natoms[0][integrator->residue_triplets[i].c]; /* psi */
                for (j = 0; j < integrator->loop_rotate_natoms[i][5]; j++)
                    integrator->loop_rotate_atoms[i][5][j] =
                        integrator->rotate_atom[0][integrator->residue_triplets[i].c][j];
            } else {
                integrator->loop_rotate_natoms[i][4] = 0;
                integrator->loop_rotate_natoms[i][5] = 0;
            }
        } else { /* psi then phi */
            j = 0;
            if (integrator->residue_triplets[i].c >= 0) {
                for (k = 0; k < top->natoms; k++)
                    if (ctx->native[k].res_num < integrator->residue_triplets[i].c) {
                        integrator
                            ->loop_rotate_atoms[i][j][integrator->loop_rotate_natoms[i][j]++] =
                            k; /* psi */
                        integrator->loop_rotate_atoms[i][j + 1]
                                                     [integrator->loop_rotate_natoms[i][j + 1]++] =
                            k; /* phi */
                    } else if (ctx->native[k].res_num == integrator->residue_triplets[i].c) {
                        if (k != ctx->native_residue[integrator->residue_triplets[i].c].O &&
                            k != ctx->native_residue[integrator->residue_triplets[i].c].C &&
                            k != ctx->native_residue[integrator->residue_triplets[i].c].CA)
                            integrator
                                ->loop_rotate_atoms[i][j][integrator->loop_rotate_natoms[i][j]++] =
                                k; /* psi */
                    }
                j += 2;
            } else {
                integrator->loop_rotate_natoms[i][4] = 0;
                integrator->loop_rotate_natoms[i][5] = 0;
            }
            if (integrator->residue_triplets[i].b >= 0) {
                for (k = 0; k < top->natoms; k++)
                    if (ctx->native[k].res_num < integrator->residue_triplets[i].b) {
                        integrator
                            ->loop_rotate_atoms[i][j][integrator->loop_rotate_natoms[i][j]++] =
                            k; /* psi */
                        integrator->loop_rotate_atoms[i][j + 1]
                                                     [integrator->loop_rotate_natoms[i][j + 1]++] =
                            k; /* phi */
                    } else if (ctx->native[k].res_num == integrator->residue_triplets[i].b) {
                        if (k != ctx->native_residue[integrator->residue_triplets[i].b].O &&
                            k != ctx->native_residue[integrator->residue_triplets[i].b].C &&
                            k != ctx->native_residue[integrator->residue_triplets[i].b].CA)
                            integrator
                                ->loop_rotate_atoms[i][j][integrator->loop_rotate_natoms[i][j]++] =
                                k; /* psi */
                    }
                j += 2;
            } else {
                integrator->loop_rotate_natoms[i][2] = 0;
                integrator->loop_rotate_natoms[i][3] = 0;
            }
            integrator->loop_rotate_natoms[i][j] =
                integrator->rotate_natoms[0][integrator->residue_triplets[i].a]; /* psi */
            for (k = 0; k < integrator->loop_rotate_natoms[i][j]; k++)
                integrator->loop_rotate_atoms[i][j][k] =
                    integrator->rotate_atom[0][integrator->residue_triplets[i].a][k];
            integrator->loop_rotate_natoms[i][j + 1] =
                integrator->rotate_natoms[1][integrator->residue_triplets[i].a]; /* phi */
            for (k = 0; k < integrator->loop_rotate_natoms[i][j + 1]; k++)
                integrator->loop_rotate_atoms[i][j + 1][k] =
                    integrator->rotate_atom[1][integrator->residue_triplets[i].a][k];
        }
    }

    for (i = 0; i < integrator->total_triplets; i++)
        for (j = 0; j < 6; j++) {
            for (k = 0; k < integrator->loop_rotate_natoms[i][j]; k++)
                integrator->loop_not_rotated[i][j][integrator->loop_rotate_atoms[i][j][k]] = 1;
            for (k = 0; k < top->natoms; k++)
                integrator->loop_not_rotated[i][j][k] = !integrator->loop_not_rotated[i][j][k];
        }

    for (i = 0; i < integrator->total_triplets; i++)
        for (j = 5; j > 0; j--) {
            if (integrator->loop_rotate_natoms[i][j] && integrator->loop_rotate_natoms[i][j - 1]) {
                for (k = 0; k < integrator->loop_rotate_natoms[i][j - 1]; k++)
                    if (integrator
                            ->loop_not_rotated[i][j][integrator->loop_rotate_atoms[i][j - 1][k]])
                        integrator->loop_int_rotate_atoms
                            [i][j - 1][integrator->loop_int_rotate_natoms[i][j - 1]++] =
                            integrator->loop_rotate_atoms[i][j - 1][k];
            }
        }

    return;
}

void InitializeSidechainRotationData(struct Topology *top, struct MCIntegrator *integrator,
                                     struct Context *ctx, struct System *sys,
                                     struct Simulation *sim) {
    int i, j, k, l, m, n;

    integrator->sidechain_torsion = (short ***)calloc(top->nresidues, sizeof(short **));
    for (i = 0; i < top->nresidues; i++) {
        integrator->sidechain_torsion[i] = (short **)calloc(
            top->amino_acids[ctx->native_residue[i].amino_num].ntorsions, sizeof(short *));
        for (j = 0; j < top->amino_acids[ctx->native_residue[i].amino_num].ntorsions; j++)
            integrator->sidechain_torsion[i][j] = (short *)calloc(4, sizeof(short));
    }

    sys->sct_E = (short *****)calloc(top->nresidues, sizeof(short ****));
    for (i = 0; i < top->nresidues; i++) {
        sys->sct_E[i] = (short ****)calloc(12, sizeof(short ***));
        for (j = 0; j < 12; j++) {
            sys->sct_E[i][j] = (short ***)calloc(12, sizeof(short **));
            for (k = 0; k < 12; k++) {
                sys->sct_E[i][j][k] = (short **)calloc(12, sizeof(short *));
                for (l = 0; l < 12; l++) {
                    sys->sct_E[i][j][k][l] = (short *)calloc(12, sizeof(short));
                }
            }
        }
    }

    for (i = 0; i < top->nresidues; i++) {
        ctx->native_residue[i].ntorsions =
            top->amino_acids[ctx->native_residue[i].amino_num].ntorsions;
        ctx->native_residue[i].nrotamers =
            top->amino_acids[ctx->native_residue[i].amino_num].nrotamers;

        /* copy rotamer angles from amino_acid structure into native_residue structure */
        for (j = 0; j < 4; j++)
            for (k = 0; k < 4; k++)
                for (l = 0; l < 4; l++)
                    for (m = 0; m < 4; m++)
                        for (n = 0; n < 4; n++)
                            ctx->native_residue[i].avg_angle[j][k][l][m][n] =
                                top->amino_acids[ctx->native_residue[i].amino_num]
                                    .avg_angle[j][k][l][m][n];

        /* rot_position gives the 4-digit base 3 representation of the rotamer */

        ctx->native_residue[i].rot_position =
            (short **)calloc(ctx->native_residue[i].nrotamers, sizeof(short *));
        for (j = 0; j < ctx->native_residue[i].nrotamers; j++)
            ctx->native_residue[i].rot_position[j] = (short *)calloc(4, sizeof(short));

        if (ctx->native_residue[i].nrotamers > 1)
            for (k = 0; k < ctx->native_residue[i].nrotamers; k++) {
                j = ctx->native_residue[i].ntorsions;
                l = k;
                do {
                    ctx->native_residue[i].rot_position[k][j - 1] = l / three[j - 1];
                    l -= three[j - 1] * ctx->native_residue[i].rot_position[k][j - 1];
                    ctx->native_residue[i].rot_position[k][j - 1] += 1;
                    j--;
                } while (j > 0);
            }
    }

    /* these structures keep track of which sidechain atoms rotate at each torsion */

    integrator->rotate_sidechain_atom = (short ***)calloc(top->nresidues, sizeof(short **));
    integrator->sidechain_not_rotated = (char ***)calloc(top->nresidues, sizeof(char **));
    for (i = 0; i < top->nresidues; i++) {
        integrator->rotate_sidechain_atom[i] =
            (short **)calloc(ctx->native_residue[i].ntorsions, sizeof(short *));
        integrator->sidechain_not_rotated[i] =
            (char **)calloc(ctx->native_residue[i].ntorsions, sizeof(char *));
        for (j = 0; j < ctx->native_residue[i].ntorsions; j++) {
            integrator->rotate_sidechain_atom[i][j] =
                (short *)calloc(top->amino_acids[ctx->native_residue[i].amino_num].rotate_natoms[j],
                                sizeof(short *));
            integrator->sidechain_not_rotated[i][j] = (char *)calloc(top->natoms, sizeof(char));
        }
    }

    integrator->rotate_sidechain_natoms = (short **)calloc(top->nresidues, sizeof(short *));
    for (i = 0; i < top->nresidues; i++)
        integrator->rotate_sidechain_natoms[i] =
            (short *)calloc(ctx->native_residue[i].ntorsions, sizeof(short));
    for (i = 0; i < top->nresidues; i++) {
        for (j = 0; j < ctx->native_residue[i].ntorsions; j++)
            integrator->rotate_sidechain_natoms[i][j] =
                top->amino_acids[ctx->native_residue[i].amino_num].rotate_natoms[j];
        for (j = 0; j < ctx->native_residue[i].ntorsions; j++)
            for (k = 0; k < 4; k++) {
                for (l = 0; l < top->natoms; l++)
                    if (ctx->native[l].res_num == i &&
                        !strcmp(ctx->native[l].atomname,
                                top->amino_acids[ctx->native_residue[i].amino_num].torsion[j][k].data())) {
                        /* sidechain_torsion records the 4 atoms that define each torsion */
                        integrator->sidechain_torsion[i][j][k] = l;
                        break;
                    }
                if (l == top->natoms) {
                    fprintf(sim->STATUS, "WARNING -- atom %s, residue %d  %s not found!\n",
                            top->amino_acids[ctx->native_residue[i].amino_num].torsion[j][k], i,
                            ctx->native_residue[i].res);
                    exit(1);
                }
            }
        for (j = 0; j < ctx->native_residue[i].ntorsions; j++)
            for (k = 0; k < integrator->rotate_sidechain_natoms[i][j]; k++) {
                for (l = 0; l < top->natoms; l++)
                    if (ctx->native[l].res_num == i &&
                        !strcmp(
                            ctx->native[l].atomname,
                            top->amino_acids[ctx->native_residue[i].amino_num].rotate_atom[j][k].data())) {
                        integrator->rotate_sidechain_atom[i][j][k] = l;
                        break;
                    }
                if (l == top->natoms) {
                    fprintf(sim->STATUS, "WARNING -- atom %s, residue %d  %s not found!\n",
                            top->amino_acids[ctx->native_residue[i].amino_num].rotate_atom[j][k], i,
                            ctx->native_residue[i].res);
                    exit(1);
                }
            }
    }

    for (i = 0; i < top->nresidues; i++)
        for (j = 0; j < ctx->native_residue[i].ntorsions; j++) {
            for (k = 0; k < top->natoms; k++)
                integrator->sidechain_not_rotated[i][j][k] = 1;
            for (k = 0; k < integrator->rotate_sidechain_natoms[i][j]; k++)
                integrator
                    ->sidechain_not_rotated[i][j][integrator->rotate_sidechain_atom[i][j][k]] = 0;
        }

    return;
}

/*=============================================================*/
/*               set up correlation data                       */
/*=============================================================*/

int SkipSelf(int s, int b, struct atom *Protein, struct residue *Residue) {
    /* returns 1 if sidechain-backbone atom pair is separated by less than 3 bonds, or 0 if not */

    if ((b == Residue[Protein[b].res_num].C || b == Residue[Protein[b].res_num].N ||
         b == Residue[Protein[b].res_num].CA) &&
        s == Residue[Protein[s].res_num].CB)
        return 1;
    else if (Residue[Protein[s].res_num].amino_num == 14)
        return 1;
    else if (b == Residue[Protein[b].res_num].CA && Protein[s].atomname[1] == 'G')
        return 1;
    else
        return 0;
}

int SkipNeighbors(int i, int j, struct atom *Protein, struct residue *Residue) {
    int first, second;

    if (Protein[i].res_num < Protein[j].res_num) {
        first  = i;
        second = j;
    } else {
        second = i;
        first  = j;
    }
    if (!(first == Residue[Protein[first].res_num].C && !strcmp(Protein[second].atomname, "CD") &&
          Residue[Protein[second].res_num].amino_num == 14) &&
        !(first == Residue[Protein[first].res_num].CA && !strcmp(Protein[second].atomname, "CD") &&
          Residue[Protein[second].res_num].amino_num == 14))
        if ((first == Residue[Protein[first].res_num].N) ||
            (second != Residue[Protein[second].res_num].CA &&
             second != Residue[Protein[second].res_num].N) ||
            Protein[first].is_sidechain || Protein[second].is_sidechain)
            return 0;

    return 1;
}

int Disulfide(int a, int b, struct atom *Protein, struct residue *Residue) {
    if ((Residue[Protein[a].res_num].amino_num == 4) &&
        (Residue[Protein[b].res_num].amino_num == 4) && !strcmp(Protein[a].atomname, "SG") &&
        !strcmp(Protein[b].atomname, "SG"))
        return 1;
    else
        return 0;
}

void CheckCorrelation(struct contact_data **Data, struct atom *Protein, struct residue *Residue,
                      int Natoms, struct Simulation *sim) {
    int i, j;

    for (i = 0; i < Natoms; i++)
        for (j = i + 1; j < Natoms; j++) {
            /* non-local residues */

            if (fabs(Protein[i].res_num - Protein[j].res_num) > sim->SKIP_LOCAL_CONTACT_RANGE) {
                if (!Disulfide(i, j, Protein, Residue)) {
                    Data[i][j].check_clashes  = 1;
                    Data[j][i].check_clashes  = 1;
                    Data[i][j].check_contacts = 1;
                    Data[j][i].check_contacts = 1;
                } else {
                    /* disulfide */
                    Data[i][j].disulfide      = 1;
                    Data[j][i].disulfide      = 1;
                    Data[i][j].check_clashes  = 0;
                    Data[j][i].check_clashes  = 0;
                    Data[i][j].check_contacts = 1;
                    Data[j][i].check_contacts = 1;
                }
            }

            /* self */

            else if (Protein[i].res_num == Protein[j].res_num) {
                if (Protein[i].is_sidechain &&
                    !Protein[j].is_sidechain) { /* sidechain - backbone */
                    if (!SkipSelf(i, j, Protein, Residue)) {
                        Data[i][j].check_clashes  = 1;
                        Data[j][i].check_clashes  = 1;
                        Data[i][j].check_contacts = 0;
                        Data[j][i].check_contacts = 0;
                    }
                } else if (!Protein[i].is_sidechain &&
                           Protein[j].is_sidechain) { /* backbone - sidechain */
                    if (!SkipSelf(j, i, Protein, Residue)) {
                        Data[i][j].check_clashes  = 1;
                        Data[j][i].check_clashes  = 1;
                        Data[i][j].check_contacts = 0;
                        Data[j][i].check_contacts = 0;
                    }
                }
            }

            /* i-i+1 */

            else if (fabs(Protein[i].res_num - Protein[j].res_num) == 1) {
                if (!SkipNeighbors(i, j, Protein, Residue)) {
                    Data[i][j].check_clashes  = 1;
                    Data[j][i].check_clashes  = 1;
                    Data[i][j].check_contacts = 0;
                    Data[j][i].check_contacts = 0;
                }
            }

            /* all other local residues */

            else if (fabs(Protein[i].res_num - Protein[j].res_num) <=
                         sim->SKIP_LOCAL_CONTACT_RANGE &&
                     fabs(Protein[i].res_num - Protein[j].res_num) >= 2) {
                Data[i][j].check_clashes  = 1;
                Data[j][i].check_clashes  = 1;
                Data[i][j].check_contacts = 0;
                Data[j][i].check_contacts = 0;
            }

            else {
                fprintf(sim->STATUS, "Uncategorizable pair\n");
                exit(0);
            }

            /* Skip relevant backbone contacts */

            if (!IsSidechainAtom(Protein[i].atomname) && !IsSidechainAtom(Protein[j].atomname) &&
                fabs(Protein[i].res_num - Protein[j].res_num) <= sim->SKIP_BB_CONTACT_RANGE) {
                Data[i][j].check_contacts = 0;
                Data[j][i].check_contacts = 0;
            }

            /* up to here, backbone contacts below the SKIP_BB_CONTACT_RANGE are turned off */
        }

    return;
}

void SetProgramOptions(struct Simulation *sim, struct System *sys, struct MCIntegrator *integrator,
                       struct Context *ctx, struct Topology *top, int argc, char *argv[]) {
    char  line[500]; /* increased to 500 from 150 */
    char  token[50];
    char  name[500]; /* increased to 500 from 50 */
    float value;
    int   find_yang_move  = 0;
    int   find_yang_scale = 0;

    std::string cfg_file;
    int  l = 0, ls = 0, MPI_STOP = 0;

    /* Acquire cfg_file name */
    cfg_file.resize(200);
    if (sim->myrank == 0) {
        strcpy(&cfg_file[0], argv[1]);
        if (argc != 2) {
            printf("ERROR!!! Usage is like this: ./fold_potential config_file, argc : %d\n", argc);
            for (l = 0; l < argc; l++) {
                printf("argc : %3d, argv : %s\n", l, argv[l]);
            }
            MPI_STOP = 1;
        }
    }

    #ifdef USE_MPI
    sim->ierr = MPI_Bcast(&MPI_STOP, 1, MPI_INT, 0, sim->mpi_world_comm);

    if (MPI_STOP == 1) {
        MPI_Finalize();
        exit(1);
    } else {
        sim->ierr = MPI_Bcast(&cfg_file[0], 200, MPI_CHAR, 0, sim->mpi_world_comm);
    }
    #endif

    /* open cfg file */
    std::ifstream file(cfg_file);
    if (!file.is_open()) {
        throw std::runtime_error("ERROR: Can't open the JSON config file: " + cfg_file);
    }

    json config = json::parse(file);
    // --- Native Protein Data ---
    std::cout << "DEBUG: Reading config file: " << cfg_file << std::endl;
    auto& npd = config["native_protein_data"];
    sim->native_file      = npd["NATIVE_FILE"].get<std::string>();
    sim->structure_file   = npd["STRUCTURE_FILE"].get<std::string>();
    sim->native_directory = npd.value("NATIVE_DIRECTORY", "None");
    sim->template_file    = npd["TEMPLATE_FILE"].get<std::string>();
    sim->alignment_file   = npd["ALIGNMENT_FILE"].get<std::string>();
    sim->pdb_out_file     = npd["PDB_OUT_FILE"].get<std::string>();
    std::cout << "DEBUG: PDB_OUT_FILE at initialization : " << sim->pdb_out_file << std::endl;
    sim->PROTEIN_NAME     = npd["PROTEIN_NAME"].get<std::string>();

    // --- Potential Parameters ---
    std::cout << "DEBUG: Reading potential parameters from config file." << std::endl;
    auto& pp = config["potential_parameters"];
    sim->NO_NEW_CLASHES       = pp["NO_NEW_CLASHES"].get<int>();
    sim->READ_POTENTIAL       = pp.value("READ_POTENTIAL", (short)0);
    sys->USE_GO_POTENTIAL     = pp["USE_GO_POTENTIAL"].get<int>();
    sys->weight_clash         = pp["CLASH_WEIGHT"].get<float>();
    sys->weight_rms           = pp["RMS_WEIGHT"].get<float>();
    // ctx->hydrogen_bond        = pp["HYDROGEN_BOND"].get<float>();
    sys->NATIVE_ATTRACTION    = pp["NATIVE_ATTRACTION"].get<float>();
    sys->NON_NATIVE_REPULSION = pp["NON_NATIVE_REPULSION"].get<float>();
    sys->NON_SPECIFIC_ENERGY  = pp["NON_SPECIFIC_ENERGY"].get<float>();

    // --- Contact Definition ---
    std::cout << "DEBUG: Reading contact definition parameters from config file." << std::endl;
    auto& cd = config["contact_definition"];
    sim->SKIP_LOCAL_CONTACT_RANGE = cd["SKIP_LOCAL_CONTACT_RANGE"].get<int>();
    sim->SKIP_BB_CONTACT_RANGE    = cd["SKIP_BB_CONTACT_RANGE"].get<int>();

    // --- Monte-Carlo Parameters ---
    sys->rmsd_constraint = config["monte_carlo_parameters"]["CONSTRAINT_RMSD"].get<float>();
    integrator->YANG_MOVE       = config["monte_carlo_parameters"]["YANG_MOVE"].get<float>();
    integrator->YANG_SCALE      = config["monte_carlo_parameters"]["YANG_SCALE"].get<int>();
    if (integrator->YANG_MOVE) { // DEBUG/NOTe: Does this turn off if YANG MOVE is not set
        find_yang_move = 1;
    }
    if (integrator->YANG_SCALE) {
        find_yang_scale = 1;
    }
    integrator->USE_SIDECHAINS = config["monte_carlo_parameters"]["USE_SIDECHAINS"].get<int>();
    sys->SEQ_DEP_HB = config["monte_carlo_parameters"]["SEQ_DEP_HB"].get<int>();
    std::cout << "DEBUG: YANG_MOVE : " << integrator->YANG_MOVE << ", YANG_SCALE : " << integrator->YANG_SCALE
              << std::endl;
    std::cout << "DEBUG: yang move flags - find_yang_move : " << find_yang_move << ", find_yang_scale : " << find_yang_scale << std::endl;

    // --- Parameter Files ---
    auto& pf = config["parameter_files"];
    sim->triplet_file          = pf["TRIPLET_ENERGY_FILE"].get<std::string>();
    sim->sctorsion_file        = pf["SIDECHAIN_TORSION_FILE"].get<std::string>();
    sim->sec_str_file          = pf["SECONDARY_STRUCTURE_FILE"].get<std::string>();
    sim->amino_data_file       = pf["AMINO_DATA_FILE"].get<std::string>();
    sim->rotamer_data_file     = pf["ROTAMER_DATA_FILE"].get<std::string>();
    sim->atom_type_file        = pf["ATOM_TYPE_FILE"].get<std::string>();
    sim->hydrogen_bonding_data = pf["HYDROGEN_BONDING_DATA"].get<std::string>();
    sim->hbond_energy_file       = pf["HYDROGEN_BOND_ENERGY_FILE"].get<std::string>();
    // sim->hbond_file            = pf["HYDROGEN_BOND_FILE"].get<std::string>();
    // sys->seq_hbond_file        = pf["SEQ_HYDROGEN_BOND_FILE"].get<std::string>();
    sim->potential_file        = pf["POTENTIAL_DATA"].get<std::string>();
    sim->aromatic_file         = pf["AROMATIC_FILE"].get<std::string>();
    // sys->all_triplet_file      = pf["TRIPLET_FILE"].get<std::string>();
    // sys->all_sctorsion_file    = pf["SC_TORISON_FILE"].get<std::string>();
    /* close cfg file */

    // --- Replica Exchange Parameters ---
    std::cout << "DEBUG: Reading replica exchange parameters from config file." << std::endl;
    auto& re = config["replica_exchange_parameters"];
    sim->NODES_PER_TEMP = re["NODES_PER_TEMP"].get<int>();
    sim->contacts_step  = re["CONTACTS_STEP"].get<int>();
    sim->umbrella       = re["UMBRELLA"].get<int>();
    integrator->MC_TEMP_MIN = re["MC_TEMP_MIN"].get<float>();
    integrator->TEMP_STEP   = re["TEMP_STEP"].get<float>();
    sys->k_bias          = re["K_BIAS"].get<float>();
    sim->number_of_contacts_max = re["NUMBER_OF_CONTACTS_MAX"].get<int>();
    top->min_seq_sep = re["MIN_SEQ_SEP"].get<int>();
    sys->contact_calpha_cutoff    = re["CONTACT_CALPHA_CUTOFF"].get<float>();

    printf("Finished reading config file!\n");
    if (strcmp(sim->native_directory.c_str(), "None") != 0) {
        char path_buf[1024];
        sprintf(path_buf, "%s/%d.pdb", sim->native_directory.c_str(), sim->myrank + sim->my_rank_offset);
        sim->native_file = path_buf;

        FILE *test_file;
        if ((test_file = fopen(sim->native_file.c_str(), "r")) == NULL) {
            int mod_rank = (sim->myrank + sim->my_rank_offset) % sim->NODES_PER_TEMP;
            // printf("The current value of mod_rank is %i", mod_rank);
            char mod_path_buf[1024];
            sprintf(mod_path_buf, "%s/%d.pdb", sim->native_directory.c_str(), mod_rank);
            sim->native_file = mod_path_buf;
            // printf("The current value of native_file is %s", native_file);
        }
    }

    std::cout << "DEBUG: number of processes : " << sim->nprocs << std::endl;
    
    if (sim->nprocs % sim->NODES_PER_TEMP != 0) {
        printf("ERROR! Number of cores must be divisible by NODES_PER_TEMP");
        //exit(1);
    }

    sim->Tnode = (float *)calloc(sim->nprocs, sizeof(float));
    sim->Enode = (float *)calloc(sim->nprocs, sizeof(float));
    sim->Cnode = (int *)calloc(sim->nprocs, sizeof(int));  // number of contacts setpoints
    sim->Nnode = (int *)calloc(sim->nprocs, sizeof(int));

    std::cout << "DEBUG: Setting up temperature and contact setpoint arrays for replica exchange." << std::endl;

    float current_temp = integrator->MC_TEMP_MIN;
    int   current_setpoint;
    for (l = 0; l < sim->nprocs; l++) {
        current_setpoint =
            sim->number_of_contacts_max - l % sim->NODES_PER_TEMP * sim->contacts_step;
        sim->Cnode[l] = current_setpoint;
        sim->Tnode[l] = current_temp;
        // current_setpoint=current_setpoint-contacts_step;
        if ((l + 1) % sim->NODES_PER_TEMP == 0) {
            current_temp = current_temp + integrator->TEMP_STEP;
        }
    }

    integrator->MC_TEMP =
        sim->Tnode[sim->myrank];  // temperature must always be externally modified if we want ot
                                  // change the sim->my_rank_offset
    sim->number_of_contacts_setpoint = sim->Cnode[sim->myrank];  /// same with setpoint

    std::cout << "DEBUG: Finished setting up temperature and contact setpoint arrays for replica exchange." << std::endl;

    /* SET name of PDB and log file  */
    std::string std_prefix = sim->pdb_out_file;
    if (sim->umbrella == 1) {
        std::string suffix = std::format("_{:5.3f}_{}", 
                                     integrator->MC_TEMP, 
                                     sim->number_of_contacts_setpoint);

        sim->std_file = std::format("{}{}.log", std_prefix, suffix);
        sim->pdb_out_file = sim->pdb_out_file + suffix;
    } else {
        sys->k_bias = 0;  // regardless of what we had inputted for k_bias, if umbrella is off, we
                          // don't want any bias!!
        if (sim->NODES_PER_TEMP == 1 && sim->my_rank_offset == 0) {
            sim->std_file = std::format("{}_T_{:5.3f}.log", 
                                    std_prefix, 
                                    integrator->MC_TEMP);

            std::string pdb_suffix = std::format("_{:5.3f}", integrator->MC_TEMP);
            sim->pdb_out_file = sim->pdb_out_file + pdb_suffix;
        } else {
            int combined_rank = sim->myrank + sim->my_rank_offset;

            sim->std_file = std::format("{}_T_{:5.3f}_{}.log", 
                                        std_prefix, 
                                        integrator->MC_TEMP, 
                                        combined_rank);

            std::string pdb_suffix = std::format("_{:5.3f}_{}", 
                                                integrator->MC_TEMP, 
                                                combined_rank);
            sim->pdb_out_file = sim->pdb_out_file + pdb_suffix;
        }
    }

    std::cout << "DEBUG: Finished setting up names for log and PDB output files." << std::endl;
    
    /* OPEN log file  for business*/
    /*sprintf(std_file, "%s_%5.3f.log", std_prefix, MC_TEMP);*/
    sim->STATUS =
        fopen(sim->std_file.c_str(), "w"); /* for some reason, code gives file handle the name STATUS */
    // fprintf(STATUS,"The value of temp1 is %s and the value of mod_rank is %i", temp1, mod_rank);
    // fflush(STATUS);

    std::cout << "DEBUG: Log file opened for writing." << std::endl;

    sim->replica_index = (int *)calloc(sim->nprocs, sizeof(int));

    // Previously:
    sim->accepted_replica  = (int *)calloc(sim->nprocs, sizeof(int));
    sim->rejected_replica  = (int *)calloc(sim->nprocs, sizeof(int));
    sim->attempted_replica = (int *)calloc(sim->nprocs, sizeof(int));

    for (l = 0; l < sim->nprocs; l++) {
        sim->accepted_replica[l] = sim->rejected_replica[l] = 0;
    }

    std::cout << "DEBUG: Replica exchange tracking arrays initialized." << std::endl;

    fprintf(sim->STATUS, "The native directory is %s \n", sim->native_directory);
    fprintf(sim->STATUS, "Temperature/setpoint range for replica exchange!\n");

    std::cout << "DEBUG: Logging initial replica exchange parameters." << std::endl;
    
    for (l = 0; l < sim->nprocs; l++) {
        fprintf(sim->STATUS, "%4d : %5.3f %d \n", l, sim->Tnode[l], sim->Cnode[l]);
    }

    fprintf(sim->STATUS, "myrank : %4d, cfg_file : %s\n k_bias : %5.3f\n setpoint : %d\n ",
            sim->myrank, cfg_file, sys->k_bias, sim->number_of_contacts_setpoint);
    fflush(sim->STATUS);

    std::cout << "DEBUG: Finished logging initial replica exchange parameters." << std::endl;

    if (find_yang_move == 0) {
        fprintf(sim->STATUS, "There is nothing on YANG_MOVE!\n");
        std::cout << "DEBUG: There is nothing on YANG_MOVE!" << std::endl;
        exit(1);
    }
    if (find_yang_scale == 0) {
        fprintf(sim->STATUS, "There is nothing on YANG_SCALE!\n");
        std::cout << "DEBUG: There is nothing on YANG_SCALE!" << std::endl;
        exit(1);
    }
    printf("Temperature/setpoint range for replica exchange!\n");

    std::cout << "DEBUG: Temperature/setpoint range for replica exchange!" << std::endl;

    /* lattice parameters */
    std::cout << "DEBUG: Setting lattice parameters." << std::endl;
    std::cout << "DEBUG: Distance dependence flag is " << DISTANCE_DEPENDENCE << std::endl;
    std::cout << "DEBUG: System lambda is " << sys->LAMBDA << " and alpha is " << sys->ALPHA << std::endl;
    if (DISTANCE_DEPENDENCE)
        integrator->LATTICE_SIZE = 1.0 / 5.25;
    else
        integrator->LATTICE_SIZE = 1 / (sys->LAMBDA * sys->ALPHA * (1.88 + 1.88));
    std::cout << "DEBUG: Lattice size set with " << sys->LAMBDA << " and " << sys->ALPHA << std::endl;
    integrator->MATRIX_SIZE      = 20;
    integrator->HALF_MATRIX_SIZE = integrator->MATRIX_SIZE / 2;
    printf("Temperature/setpoint range for replica exchange2!\n");

    std::cout << "DEBUG: Finished setting lattice parameters." << std::endl;

    /*Read the constraints--AB*/
    if (sim->constraint_file[0] != '\0' && strcmp(sim->constraint_file.c_str(), "None") != 0) {
        fprintf(sim->STATUS, "There are constraints to read!\n");
        fflush(sim->STATUS);
        Read_constraints(sys, sim, ctx);
        fprintf(sim->STATUS, "Constraints read\n");
        fflush(sim->STATUS);
    } else {
        // Optional: Log that you are skipping constraints
        fprintf(sim->STATUS, "No constraints provided. Skipping.\n");
        fflush(sim->STATUS);
    }
    printf("Finished setting program options!\n");

    std::cout << "DEBUG: Finished setting program options." << std::endl;

    return;
}