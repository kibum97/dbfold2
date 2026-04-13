#include "energy.h"
#include "define.h"
#include "constraint.h"
#include "contacts.h"
#include "init.h"
#include "loop.h"
#include "vector.h"
#include "hbonds.h"
#include <string_view>

void ResetEnergies(
    struct Context *ctx,
    struct Simulation *sim,
    struct System *sys,
    struct Topology *top,
    long int check) {
    TypeContacts(ctx, sys, top);
    ctx->E_pot   = FullAtomEnergy(sys, ctx);
    ctx->E_tor   = torsionenergy(ctx, top, sys);
    ctx->E_sct   = sctenergy(ctx, sys, top);
    ctx->E_aro   = aromaticenergy(ctx, sys, top);
    ctx->E_hbond = HydrogenBonds(top, ctx, sys);

    if (sim->constraint_file != "None") {                            // AB
        ctx->E_constraint = Compute_constraint_energy(ctx, sys, ctx->native_residue, ctx->native);  // AB
    } else if (sim->rmsd_constraint_file != "None") {
        fprintf(sim->STATUS, "RMSD constraint file is not None \n");
        ctx->E_constraint =
            Compute_frag_rmsd_constraint_energy(sys, ctx->struct_f1, ctx->struct_f2, constraint_align);  // AB
    } else {                                                                              // AB
        ctx->E_constraint = 0;                                                                 // AB
    }  // AB

    ctx->E = sys->weight_potential * ctx->E_pot + sys->weight_clash * ctx->nclashes + sys->weight_hbond * ctx->E_hbond +
        TOR_WEIGHT * ctx->E_tor + SCT_WEIGHT * ctx->E_sct + ARO_WEIGHT * ctx->E_aro + ctx->E_constraint;

    if (abs(ctx->E - ctx->prev_E) > 1.0 && check) {
        fprintf(sim->STATUS, "ResetEnergies(%ld): prev_E(%.5f), New E(%.5f) dE = %.5f\n", check, ctx->prev_E,
                ctx->E, ctx->E - ctx->prev_E);
        fprintf(sim->STATUS,
                "E_pot %.5f %.5f, E_hbond %.5f %.5f, E_tor %.5f %.5f, E_sct %.5f %.5f, E_aro %.5f "
                "%.5f E_constraint %.5f %.5f \n",
                ctx->prev_E_pot, ctx->E_pot, ctx->prev_E_hbond, ctx->E_hbond, ctx->prev_E_tor, ctx->E_tor, ctx->prev_E_sct, ctx->E_sct,
                ctx->prev_E_aro, ctx->E_aro, ctx->prev_E_constraint, ctx->E_constraint);
    }

    ctx->prev_E_pot   = ctx->E_pot;
    ctx->prev_E_tor   = ctx->E_tor;
    ctx->prev_E_sct   = ctx->E_sct;
    ctx->prev_E_aro   = ctx->E_aro;
    ctx->prev_E_hbond = ctx->E_hbond;

    ctx->prev_E_constraint = ctx->E_constraint;  // AB

    ctx->prev_E = ctx->E;
}

float FullAtomEnergy(
    struct System *sys,
    struct Context *ctx
) {
    float e = 0;
    int   i, j;

    for (i = 0; i < sys->MAX_TYPES; i++) {
        for (j = i; j < sys->MAX_TYPES; j++) {
            e += ctx->type_contacts[i][j] * sys->potential[i][j];
            //      if(((i==bb_O_type)&&(j==bb_N_type))||((i==bb_N_type)&&(j==bb_O_type)))
            //	fprintf(STATUS, "%d %d %d %5.2f\n", i, j, type_contacts[i][j], potential[i][j]);
        }
    }
    return e;
}

float torsionenergy(
    struct Context *ctx,
    struct Topology *top,
    struct System *sys
) {
    int   cur_res, x, y, z, w;
    float locang_int = 30.;
    float loctor_int = 60.;

    float e = 0;
    int   i;

    check_bb(top, ctx);
    for (i = 0; i < top->nresidues - 2; i++) {
        cur_res = i;
        x       = (int)((ctx->a_PCA[i]) / locang_int);
        y       = (int)((ctx->a_bCA[i]) / locang_int);
        z       = (int)((ctx->phim[i]) / loctor_int);
        w       = (int)((ctx->psim[i]) / loctor_int);
        e += sys->torsion_E[cur_res][x][y][z][w];
    }

    e /= 1000.;
    return e;
}

float sctenergy(
    struct Context *ctx,
    struct System *sys,
    struct Topology *top
) {
    int   cur_res;
    int   i_ang[4];
    float locang_int = 30.;
    float ang[4];

    float e = 0;
    int   i, j;

    for (i = 0; i < top->nresidues - 2; i++) {
        if (ctx->native_residue[i + 1].ntorsions == 0)
            continue;
        cur_res = i;
        for (j = 0; j < 4; j++)
            i_ang[j] = 0;
        for (j = 0; j < ctx->native_residue[i + 1].ntorsions; j++) {
            ang[j] = ctx->native_residue[i + 1].tmpchi[j] * rad2deg;
            if (ang[j] < -180.)
                while (ang[j] < -180.)
                    ang[j] += 360.;
            else if (ang[j] >= 180.)
                while (ang[j] >= 180.)
                    ang[j] -= 360.;
            if (ang[j] < 0)
                ang[j] += +180. + 0.00001;
            else
                ang[j] += +180. - 0.00001;
            ang[j] += 15.;
            ang[j] = fmod(ang[j], 360.);
            //        fprintf(STATUS, "%3d %3d %9.3f\n", i, j, ang[j]);
            i_ang[j] = (int)(ang[j] / locang_int);
        }
        e += sys->sct_E[cur_res][i_ang[0]][i_ang[1]][i_ang[2]][i_ang[3]];
    }

    e /= 1000.;
    return e;
}

float aromaticenergy(
    struct Context *ctx,
    struct System *sys,
    struct Topology *top
) {
    float         e = 0;
    int           r1, r2;
    int           i, j, x;
    Vec3 vector_1, vector_2, vector_center;
    Vec3 plane_1, plane_2;
    float         distance_2, plane_angle;
    float         aro_int = 10.;

    //  if(mcstep>=1)
    //    exit(0);
    for (r1 = 0; r1 < top->Naromatic; r1++)
        for (r2 = r1 + 1; r2 < top->Naromatic; r2++) {
            i = top->Res_aromatic[r1];
            j = top->Res_aromatic[r2];
            aromatic_center(ctx, i, &vector_1);
            aromatic_center(ctx, j, &vector_2);
            MakeVector(vector_1, vector_2, &vector_center);
            distance_2 = D2(vector_1, vector_2);
            if (distance_2 < (AROMATIC_DISTANCE * AROMATIC_DISTANCE)) {
                plane_angle = aromatic_plane(ctx, i, j, &plane_1, &plane_2);
                if (plane_angle > 89.9)
                    plane_angle = 89.9;
                x = (int)(plane_angle / aro_int);
                e += sys->aromatic_E[x];
                //        fprintf(STATUS, "%3d %3d %3d %3d %8.3f %2d %5d %8.3f\n", r1, r2, i, j,
                //        plane_angle, x, aromatic_E[x], e);
            }
        }

    e /= 1000.;
    return e;
}

void aromatic_center(
    struct Context *ctx,
    int res_no, Vec3 *V) {
    const auto& res = ctx->native_residue[res_no];
    const auto& atoms = ctx->native;
    
    if ((strcmp(res.res, "PHE") == 0) || (strcmp(res.res, "RING") == 0)) {
        *V = (atoms[res.CG].xyz + atoms[res.CE1].xyz + atoms[res.CE2].xyz) / 3.;
    } else if (strcmp(res.res, "TRP") == 0) {
        *V = (atoms[res.CG].xyz + atoms[res.CZ2].xyz + atoms[res.CZ3].xyz) / 3.;
    }
    return;
}

float aromatic_plane(
    struct Context *ctx,
    int res_a, int res_b, Vec3 *plane_a, Vec3 *plane_b) {

    const auto& atoms = ctx->native;
    const auto& res_data_a = ctx->native_residue[res_a];
    const auto& res_data_b = ctx->native_residue[res_b];
    
    std::string_view name_a = res_data_a.res;
    std::string_view name_b = res_data_b.res;

    Vec3 v1, v2;
    if (name_a == "PHE" || name_a == "RING") {
        v1 = atoms[res_data_a.CE1].xyz - atoms[res_data_a.CG].xyz;
        v2 = atoms[res_data_a.CE2].xyz - atoms[res_data_a.CG].xyz;
    } else if (name_a == "TRP") {
        v1 = atoms[res_data_a.CZ2].xyz - atoms[res_data_a.CG].xyz;
        v2 = atoms[res_data_a.CZ3].xyz - atoms[res_data_a.CG].xyz;
    }
    *plane_a = v1.cross(v2); // Direct Eigen cross product

    Vec3 v3, v4;
    if (name_b == "PHE" || name_b == "RING") {
        v3 = atoms[res_data_b.CE1].xyz - atoms[res_data_b.CG].xyz;
        v4 = atoms[res_data_b.CE2].xyz - atoms[res_data_b.CG].xyz;
    } else if (name_b == "TRP") {
        v3 = atoms[res_data_b.CZ2].xyz - atoms[res_data_b.CG].xyz;
        v4 = atoms[res_data_b.CZ3].xyz - atoms[res_data_b.CG].xyz;
    }
    *plane_b = v3.cross(v4);

    // Angle calculation
    float angle = Angle(*plane_a, *plane_b) * rad2deg;
    if (angle > 90.0f) 
        angle = 180.0f - angle;

    return angle;
}
