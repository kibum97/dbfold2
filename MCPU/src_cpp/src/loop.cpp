#include "loop.h"

#include <stdio.h>
#include <span>
#include <array>
#include <fstream>
#include <sstream>
#include <string>

#include "contacts.h"
#include "define.h"

#include "init.h"
#include "jac_local.h"
#include "lattice_util.h"
#include "misc_util.h"
#include "tripep_closure.h"
#include "vector.h"

//========================================================================================================
void Read_movable_region(
    struct Simulation *sim,
    struct Topology *top
) {
    // This function reads in the movable region from a file
    char  str[1000];
    int   movable_lower_bound;
    int   movable_upper_bound;
    int   n_movable = 0;
    FILE *fp        = fopen(sim->movable_region_file.c_str(), "r");
    if (fp == NULL) {
        fprintf(sim->STATUS, "Could not open file %s \n", sim->movable_region_file);
        exit(1);
    } else {
        fprintf(sim->STATUS, "Successfully opened movable region file %s \n", sim->movable_region_file);
        fflush(sim->STATUS);
    }
    while (fgets(str, 1000, fp) != NULL) {
        if (str[0] != '#') {
            if (sscanf(str, "%d %d", &movable_lower_bound, &movable_upper_bound) != 2)
                continue;
            for (int i = movable_lower_bound; i < movable_upper_bound; i++) {
                top->movable_residue_map[n_movable] = i;
                n_movable++;
            }
        }
    }
    fclose(fp);
    for (int i = 0; i < n_movable; i++) {
        fprintf(sim->STATUS, "DEBUG:Movable residue %d is %d \n", i, top->movable_residue_map[i]);
        fflush(sim->STATUS);
    }
    top->n_movable_residues = n_movable;
    return;
}

void integloop(
    struct Context *ctx, struct Topology *top, struct MCIntegrator *integrator, struct System *sys, struct Simulation *sim,
    float step_size, int *n_soln)  // Does a local move...based on the Dill paper
{
    Eigen::Matrix<double, 3, 5> r_n, r_a, r_c, r_o;
    std::array<Eigen::Matrix<double, 3, MSAR>, 5> r_s;
    Eigen::Matrix<double, 3, 5> rorg_n, rorg_a, rorg_c, rorg_o;
    std::array<Eigen::Matrix<double, 3, MSAR>, 5> rorg_s;
    char   out_pdb[200];
    char   res_name[5][4];
    int    write_out_pdb = 0;

    double     b_len[6], b_ang[7];
    double     t_len[6], t_ang[7];
    double     o_len[6], o_ang[7];
    double     small_value = 0.001;
    int        ns[5];
    int        sel_res, n0;
    static int step;
    static int success = 0;
    //--------------------------------------------------------------------------------------------------------
    *n_soln = 0;  // Number of solutions to 16th degree polynomial
    step++;

    ctx->mc.loop_size       = 1;
    ctx->all_rotated_natoms = 0;
    ctx->total_pairs = ctx->total_pairs2 = 0;
    ctx->total_hbond_pairs          = 0;

    top->res_atomno[top->nresidues] = top->natoms;

    sel_res = (int)(threefryrand() * top->n_movable_residues);  // Select a residue to initialize local move
    sel_res = top->movable_residue_map[sel_res];          // KP editted for movable region constraint
                                                     //  sel_res = 20;
    ctx->mc.is_phi = (int)(threefryrand() * 2);  // With 50% probability, rotate about phi angle, otherwise psi angle
    step_size = integrator->YANG_SCALE * 2. * deg2rad * GaussianNum();  // Magnitude  of driver rotation, in radians
    

    if (ctx->mc.is_phi)
    // residues are to be formed along (-) direction.
    {
        n0 = sel_res - 3;
        if (sel_res <
            4)  // Don't make a move, since you don't have enough residues before to work with
        {
            //      fprintf(sim->STATUS, "Too left! n0: %d\n", n0);
            return;
        }
        if ((ctx->native_residue[sel_res].amino_num == 14) ||
            (ctx->native_residue[sel_res - 1].amino_num == 14) ||
            (ctx->native_residue[sel_res - 2].amino_num == 14) ||
            (ctx->native_residue[sel_res - 3].amino_num == 14)) {
            //      fprintf(sim->STATUS, "Can't rotate phi angle of proline! n0: %d\n", n0);
            return;
        }
        if ((top->is_template[sel_res] == 1) || (top->is_template[sel_res - 1] == 1) ||
            (top->is_template[sel_res - 2] == 1) || (top->is_template[sel_res - 3] == 1))
            return;
    } else
    // residues are to be formed along (+) direction.
    {
        n0 = sel_res + 1;
        if (sel_res > top->nresidues - 5) {
            //      fprintf(STATUS, "Too right! n0: %d\n", n0);
            return;
        }
        if ((ctx->native_residue[sel_res + 1].amino_num == 14) ||
            (ctx->native_residue[sel_res + 2].amino_num == 14) ||
            (ctx->native_residue[sel_res + 3].amino_num == 14)) {
            //      fprintf(STATUS, "Can't rotate phi angle of proline! n0: %d\n", n0);
            return;
        }
        if ((top->is_template[sel_res] == 1) || (top->is_template[sel_res + 1] == 1) ||
            (top->is_template[sel_res + 2] == 1) || (top->is_template[sel_res + 3] == 1))
            return;
    }
    //  fprintf(STATUS, "Local move at %d\n", n0);
    strcpy(res_name[0], ctx->native_residue[n0 - 1].res);
    strcpy(res_name[1], ctx->native_residue[n0].res);
    strcpy(res_name[2], ctx->native_residue[n0 + 1].res);
    strcpy(res_name[3], ctx->native_residue[n0 + 2].res);
    strcpy(res_name[4], ctx->native_residue[n0 + 3].res);

    get_orgco(top, ctx, rorg_n, rorg_a, rorg_c, rorg_o, rorg_s, ns, n0);
    get_coord(top, ctx, r_n, r_a, r_c, r_o, r_s, ns, n0);

    c_bnd_len(rorg_a.col(1), rorg_c.col(1), &o_len[0]);
    c_bnd_len(rorg_c.col(1), rorg_n.col(2), &o_len[1]);
    c_bnd_len(rorg_n.col(2), rorg_a.col(2), &o_len[2]);
    c_bnd_len(rorg_a.col(2), rorg_c.col(2), &o_len[3]);
    c_bnd_len(rorg_c.col(2), rorg_n.col(3), &o_len[4]);
    c_bnd_len(rorg_n.col(3), rorg_a.col(3), &o_len[5]);
    c_bnd_ang(rorg_n.col(1), rorg_a.col(1), rorg_c.col(1), &o_ang[0]);
    c_bnd_ang(rorg_a.col(1), rorg_c.col(1), rorg_n.col(2), &o_ang[1]);
    c_bnd_ang(rorg_c.col(1), rorg_n.col(2), rorg_a.col(2), &o_ang[2]);
    c_bnd_ang(rorg_n.col(2), rorg_a.col(2), rorg_c.col(2), &o_ang[3]);
    c_bnd_ang(rorg_a.col(2), rorg_c.col(2), rorg_n.col(3), &o_ang[4]);
    c_bnd_ang(rorg_c.col(2), rorg_n.col(3), rorg_a.col(3), &o_ang[5]);
    c_bnd_ang(rorg_n.col(3), rorg_a.col(3), rorg_c.col(3), &o_ang[6]);

    c_bnd_len(r_a.col(1), r_c.col(1), &b_len[0]);
    c_bnd_len(r_c.col(1), r_n.col(2), &b_len[1]);
    c_bnd_len(r_n.col(2), r_a.col(2), &b_len[2]);
    c_bnd_len(r_a.col(2), r_c.col(2), &b_len[3]);
    c_bnd_len(r_c.col(2), r_n.col(3), &b_len[4]);
    c_bnd_len(r_n.col(3), r_a.col(3), &b_len[5]);
    c_bnd_ang(r_n.col(1), r_a.col(1), r_c.col(1), &b_ang[0]);
    c_bnd_ang(r_a.col(1), r_c.col(1), r_n.col(2), &b_ang[1]);
    c_bnd_ang(r_c.col(1), r_n.col(2), r_a.col(2), &b_ang[2]);
    c_bnd_ang(r_n.col(2), r_a.col(2), r_c.col(2), &b_ang[3]);
    c_bnd_ang(r_a.col(2), r_c.col(2), r_n.col(3), &b_ang[4]);
    c_bnd_ang(r_c.col(2), r_n.col(3), r_a.col(3), &b_ang[5]);
    c_bnd_ang(r_n.col(3), r_a.col(3), r_c.col(3), &b_ang[6]);

    // Now we call yangloop, which actually does the rotation

    yangloop(ctx, sim, r_n, r_a, r_c, r_o, r_s, ns, b_len, b_ang, n0, step_size, res_name, n_soln);
    c_bnd_len(r_a.col(1), r_c.col(1), &t_len[0]);
    c_bnd_len(r_c.col(1), r_n.col(2), &t_len[1]);
    c_bnd_len(r_n.col(2), r_a.col(2), &t_len[2]);
    c_bnd_len(r_a.col(2), r_c.col(2), &t_len[3]);
    c_bnd_len(r_c.col(2), r_n.col(3), &t_len[4]);
    c_bnd_len(r_n.col(3), r_a.col(3), &t_len[5]);
    c_bnd_ang(r_n.col(1), r_a.col(1), r_c.col(1), &t_ang[0]);
    c_bnd_ang(r_a.col(1), r_c.col(1), r_n.col(2), &t_ang[1]);
    c_bnd_ang(r_c.col(1), r_n.col(2), r_a.col(2), &t_ang[2]);
    c_bnd_ang(r_n.col(2), r_a.col(2), r_c.col(2), &t_ang[3]);
    c_bnd_ang(r_a.col(2), r_c.col(2), r_n.col(3), &t_ang[4]);
    c_bnd_ang(r_c.col(2), r_n.col(3), r_a.col(3), &t_ang[5]);
    c_bnd_ang(r_n.col(3), r_a.col(3), r_c.col(3), &t_ang[6]);
    if (*n_soln > 0) {
        if (fabs(t_len[0] - o_len[0]) > small_value) {
            //      fprintf(STATUS, "distance 0: %6.3lf %6.3lf\n", t_len[0], o_len[0]);
            return;
        }
        if (fabs(t_len[1] - o_len[1]) > small_value) {
            //      fprintf(STATUS, "distance 1: %6.3lf %6.3lf\n", t_len[1], o_len[1]);
            return;
        }
        if (fabs(t_len[2] - o_len[2]) > small_value) {
            //      fprintf(STATUS, "distance 2: %6.3lf %6.3lf\n", t_len[2], o_len[2]);
            return;
        }
        if (fabs(t_len[3] - o_len[3]) > small_value) {
            //      fprintf(STATUS, "distance 3: %6.3lf %6.3lf\n", t_len[3], o_len[3]);
            return;
        }
        if (fabs(t_len[4] - o_len[4]) > small_value) {
            //      fprintf(STATUS, "distance 4: %6.3lf %6.3lf\n", t_len[4], o_len[4]);
            return;
        }
        if (fabs(t_len[5] - o_len[5]) > small_value) {
            //      fprintf(STATUS, "distance 5: %6.3lf %6.3lf\n", t_len[5], o_len[5]);
            return;
        }

        if (fabs(t_ang[0] - o_ang[0]) > small_value * deg2rad * 100.) {
            //      fprintf(STATUS, "angle 0:  %6.3lf %6.3lf\n", t_ang[0]*rad2deg,
            //      o_ang[0]*rad2deg);
            return;
        }
        if (fabs(t_ang[1] - o_ang[1]) > small_value * deg2rad * 100.) {
            //      fprintf(STATUS, "angle 1:  %6.3lf %6.3lf\n", t_ang[1]*rad2deg,
            //      o_ang[1]*rad2deg);
            return;
        }
        if (fabs(t_ang[2] - o_ang[2]) > small_value * deg2rad * 100.) {
            //      fprintf(STATUS, "angle 2:  %6.3lf %6.3lf\n", t_ang[2]*rad2deg,
            //      o_ang[2]*rad2deg);
            return;
        }
        if (fabs(t_ang[3] - o_ang[3]) > small_value * deg2rad * 100.) {
            //      fprintf(STATUS, "angle 3:  %6.3lf %6.3lf\n", t_ang[3]*rad2deg,
            //      o_ang[3]*rad2deg);
            return;
        }
        if (fabs(t_ang[4] - o_ang[4]) > small_value * deg2rad * 100.) {
            //      fprintf(STATUS, "angle 4:  %6.3lf %6.3lf\n", t_ang[4]*rad2deg,
            //      o_ang[4]*rad2deg);
            return;
        }
        if (fabs(t_ang[5] - o_ang[5]) > small_value * deg2rad * 100.) {
            //      fprintf(STATUS, "angle 5:  %6.3lf %6.3lf\n", t_ang[5]*rad2deg,
            //      o_ang[5]*rad2deg);
            return;
        }
        if (fabs(t_ang[6] - o_ang[6]) > small_value * deg2rad * 100.) {
            //      fprintf(STATUS, "angle 6:  %6.3lf %6.3lf\n", t_ang[6]*rad2deg,
            //      o_ang[6]*rad2deg);
            return;
        }
        //    fprintf(STATUS, "passed\n");
    }

    if (*n_soln > 0) {
        put_coord(top, ctx, r_n, r_a, r_c, r_o, r_s, ns, n0);
        yang_rotate(top, integrator, ctx, n0, ctx->mc.is_phi);
        ctx->all_rotated_natoms = integrator->yang_rotated_natoms;
        ctx->all_rotated_atoms  = integrator->yang_rotated_atoms;
        UpdateLattice(ctx, integrator, ctx->all_rotated_natoms, ctx->all_rotated_atoms);
        NewDeltaContacts(ctx, sys, ctx->all_rotated_natoms, ctx->all_rotated_atoms, integrator->yang_not_rotated);
        success++;
    }

    if (write_out_pdb) {
        sprintf(out_pdb, "data/135l_loop_min.pdb");
        fprintf(sim->STATUS, "\nRecording the solution at min. rmsd in %s\n", out_pdb);
        write_pdb_backbone(sim, out_pdb, res_name, r_n, r_a, r_c, r_o, r_s, n0, n0 + 4);
    }
    // fprintf(STATUS,  "Yang move attempted \n");
    return;
}

//========================================================================================================
void get_orgco(
    struct Topology *top, struct Context *ctx,
    Eigen::Matrix<double, 3, 5>& rorg_n,
    Eigen::Matrix<double, 3, 5>& rorg_a,
    Eigen::Matrix<double, 3, 5>& rorg_c,
    Eigen::Matrix<double, 3, 5>& rorg_o,
    std::array<Eigen::Matrix<double, 3, MSAR>, 5>& rorg_s,
    int ns[5], int n0) {

    for (int i = 0; i < 5; i++)
        ns[i] = 0;

    for (int cur_res = 0; cur_res < 5; cur_res++)
        for (int cur_atom = top->res_atomno[n0 - 1]; cur_atom < top->res_atomno[n0 + 4]; cur_atom++)
            if (ctx->native[cur_atom].res_num == n0 - 1 + cur_res)
                get_orgat(ctx, rorg_n, rorg_a, rorg_c, rorg_o, rorg_s, ns, cur_res, cur_atom);

    return;
}

//========================================================================================================
void get_orgat(
    struct Context *ctx,
    Eigen::Matrix<double, 3, 5>& rorg_n,   // Pass by reference (&)
    Eigen::Matrix<double, 3, 5>& rorg_a,
    Eigen::Matrix<double, 3, 5>& rorg_c,
    Eigen::Matrix<double, 3, 5>& rorg_o,
    std::array<Eigen::Matrix<double, 3, MSAR>, 5>& rorg_s, 
    int ns[5], int r, int a){

    std::string_view name = ctx->native[a].atomname;
    const Vec3& pos = ctx->orig_native[a].xyz;

    if (name == "N") {
        rorg_n.col(r) = pos;
    } else if (name == "CA") {
        rorg_a.col(r) = pos;
    } else if (name == "C") {
        rorg_c.col(r) = pos;
    } else if (name == "O") {
        rorg_o.col(r) = pos;
    } else {
        // Assign to the sidechain matrix for residue 'r' at column 'ns[r]'
        rorg_s[r].col(ns[r]) = pos;
        ns[r]++;
    }
}
//========================================================================================================
void get_coord(
    struct Topology *top, struct Context *ctx,
    Eigen::Matrix<double, 3, 5>& r_n,
    Eigen::Matrix<double, 3, 5>& r_a,
    Eigen::Matrix<double, 3, 5>& r_c,
    Eigen::Matrix<double, 3, 5>& r_o,
    std::array<Eigen::Matrix<double, 3, MSAR>, 5>& r_s,
    int ns[5], int n0) {
    for (int i = 0; i < 5; i++)
        ns[i] = 0;

    for (int cur_res = 0; cur_res < 5; cur_res++)
        for (int cur_atom = top->res_atomno[n0 - 1]; cur_atom < top->res_atomno[n0 + 4]; cur_atom++)
            if (ctx->native[cur_atom].res_num == n0 - 1 + cur_res)
                get_atom(ctx, r_n, r_a, r_c, r_o, r_s, ns, cur_res, cur_atom);

    return;
}

//========================================================================================================
void get_atom(
    struct Context *ctx,
    Eigen::Matrix<double, 3, 5>& r_n,
    Eigen::Matrix<double, 3, 5>& r_a,
    Eigen::Matrix<double, 3, 5>& r_c,
    Eigen::Matrix<double, 3, 5>& r_o,
    std::array<Eigen::Matrix<double, 3, MSAR>, 5>& r_s,
    int ns[5], int r, int a) {

    std::string_view name = ctx->native[a].atomname;
    const Vec3& pos = ctx->orig_native[a].xyz;
    
    if (name == "N") {
        r_n.col(r) = pos;
    } else if (name == "CA") {
        r_a.col(r) = pos;
    } else if (name == "C") {
        r_c.col(r) = pos;
    } else if (name == "O") {
        r_o.col(r) = pos;
    } else {
        r_s[r].col(ns[r]) = pos;
        ns[r]++;
    }
    return;
}

//========================================================================================================
void put_coord(
    struct Topology *top, struct Context *ctx,
    const Eigen::Matrix<double, 3, 5>& r_n,
    const Eigen::Matrix<double, 3, 5>& r_a,
    const Eigen::Matrix<double, 3, 5>& r_c,
    const Eigen::Matrix<double, 3, 5>& r_o,
    const std::array<Eigen::Matrix<double, 3, MSAR>, 5>& r_s,
    int ns[5], int n0) {

    for (int i = 0; i < 5; i++)
        ns[i] = 0;

    for (int cur_res = 0; cur_res < 5; cur_res++)
        for (int cur_atom = top->res_atomno[n0 - 1]; cur_atom < top->res_atomno[n0 + 4]; cur_atom++)
            if (ctx->native[cur_atom].res_num == n0 - 1 + cur_res)
                put_atom(ctx, r_n, r_a, r_c, r_o, r_s, ns, cur_res, cur_atom);
    return;
}

//========================================================================================================
void put_atom(
    struct Context *ctx,
    const Eigen::Matrix<double, 3, 5>& r_n,
    const Eigen::Matrix<double, 3, 5>& r_a,
    const Eigen::Matrix<double, 3, 5>& r_c,
    const Eigen::Matrix<double, 3, 5>& r_o,
    const std::array<Eigen::Matrix<double, 3, MSAR>, 5>& r_s,
    int ns[5], int r, int a) {
    
    std::string_view name = ctx->native[a].atomname;
    Vec3& pos = ctx->orig_native[a].xyz;
    
    if (name == "N") {
        pos = r_n.col(r);
    } else if (name == "CA") {
        pos = r_a.col(r);
    } else if (name == "C") {
        pos = r_c.col(r);
    } else if (name == "O") {
        pos = r_o.col(r);
    } else {
        pos = r_s[r].col(ns[r]);
        ns[r]++;
    }
    return;
}
//========================================================================================================
void yang_rotate(
    struct Topology *top, struct MCIntegrator *integrator, struct Context *ctx,
    int n0, int is_phi) {
    int i, cur_atom;
    int another = 0;

    for (i = 0; i < top->natoms; i++)
        integrator->yang_not_rotated[i] = 1;

    is_rotated[ctx->native_residue[n0].C]            = 1;
    integrator->yang_not_rotated[ctx->native_residue[n0].C]      = 0;
    integrator->yang_rotated_atoms[another++]               = ctx->native_residue[n0].C;
    is_rotated[ctx->native_residue[n0].O]            = 2;
    integrator->yang_not_rotated[ctx->native_residue[n0].O]      = 0;
    integrator->yang_rotated_atoms[another++]               = ctx->native_residue[n0].O;
    is_rotated[ctx->native_residue[n0 + 1].N]        = 2;
    integrator->yang_not_rotated[ctx->native_residue[n0 + 1].N]  = 0;
    integrator->yang_rotated_atoms[another++]               = ctx->native_residue[n0 + 1].N;
    is_rotated[ctx->native_residue[n0 + 1].CA]       = 3;
    integrator->yang_not_rotated[ctx->native_residue[n0 + 1].CA] = 0;
    integrator->yang_rotated_atoms[another++]               = ctx->native_residue[n0 + 1].CA;
    is_rotated[ctx->native_residue[n0 + 1].C]        = 4;
    integrator->yang_not_rotated[ctx->native_residue[n0 + 1].C]  = 0;
    integrator->yang_rotated_atoms[another++]               = ctx->native_residue[n0 + 1].C;
    is_rotated[ctx->native_residue[n0 + 1].O]        = 4;
    integrator->yang_not_rotated[ctx->native_residue[n0 + 1].O]  = 0;
    integrator->yang_rotated_atoms[another++]               = ctx->native_residue[n0 + 1].O;
    is_rotated[ctx->native_residue[n0 + 2].N]        = 5;
    integrator->yang_not_rotated[ctx->native_residue[n0 + 2].N]  = 0;
    integrator->yang_rotated_atoms[another++]               = ctx->native_residue[n0 + 2].N;

    if (is_phi) {
        is_rotated[ctx->native_residue[n0 + 2].CA]       = 6;
        integrator->yang_not_rotated[ctx->native_residue[n0 + 2].CA] = 0;
        integrator->yang_rotated_atoms[another++]               = ctx->native_residue[n0 + 2].CA;
        is_rotated[ctx->native_residue[n0 + 2].C]        = 6;
        integrator->yang_not_rotated[ctx->native_residue[n0 + 2].C]  = 0;
        integrator->yang_rotated_atoms[another++]               = ctx->native_residue[n0 + 2].C;
        is_rotated[ctx->native_residue[n0 + 2].O]        = 6;
        integrator->yang_not_rotated[ctx->native_residue[n0 + 2].O]  = 0;
        integrator->yang_rotated_atoms[another++]               = ctx->native_residue[n0 + 2].O;
    } else {
        is_rotated[ctx->native_residue[n0 - 1].O]       = 6;
        integrator->yang_not_rotated[ctx->native_residue[n0 - 1].O] = 0;
        integrator->yang_rotated_atoms[another++]              = ctx->native_residue[n0 - 1].O;
        is_rotated[ctx->native_residue[n0].N]           = 6;
        integrator->yang_not_rotated[ctx->native_residue[n0].N]     = 0;
        integrator->yang_rotated_atoms[another++]              = ctx->native_residue[n0].N;
        is_rotated[ctx->native_residue[n0].CA]          = 6;
        integrator->yang_not_rotated[ctx->native_residue[n0].CA]    = 0;
        integrator->yang_rotated_atoms[another++]              = ctx->native_residue[n0].CA;
    }

    for (cur_atom = top->res_atomno[n0]; cur_atom < top->res_atomno[n0 + 3]; cur_atom++) {
        if (ctx->native[cur_atom].is_sidechain) {
            is_rotated[cur_atom]          = 2 * (ctx->native[cur_atom].res_num - n0) + 1;
            integrator->yang_not_rotated[cur_atom]    = 0;
            integrator->yang_rotated_atoms[another++] = cur_atom;
        }
    }

    integrator->yang_rotated_natoms = another;

    return;
}
//========================================================================================================
void yangloop(
    struct Context *ctx, struct Simulation *sim,
    Eigen::Matrix<double, 3, 5>& r_n,
    Eigen::Matrix<double, 3, 5>& r_a,
    Eigen::Matrix<double, 3, 5>& r_c,
    Eigen::Matrix<double, 3, 5>& r_o,
    std::array<Eigen::Matrix<double, 3, MSAR>, 5>& r_s,
    int ns[5], double b_len[6], double b_ang[7], int n0,
    float dih_ch, char res_name[5][4], int *n_soln) {
    
    int    i, j, k, z;
    double rmsd, sum, dr[3];
    int    calc_rmsd = 1;
    double t_ang[2];
    int    index         = 1000;
    double min_value     = 1000.;
    int    n_soln_before = 0;
    double before_jacobi, after_jacobi;
    
    std::array<Eigen::Matrix3d, max_soln> r_soln_n;
    std::array<Eigen::Matrix3d, max_soln> r_soln_a;
    std::array<Eigen::Matrix3d, max_soln> r_soln_c;
    std::array<Eigen::Matrix3d, max_soln> r_soln_o;
    std::array<Eigen::Matrix<double, 3, MSAR>, 3> r_soln_s;

    
    Eigen::Matrix<double, 3, 3> r0_n, r0_a, r0_c, r0_o;
    // std::array<Eigen::Matrix<double, 3, MSAR>, 3> r0_s;
    Eigen::Matrix<double, 3, 3> r0drms_n, r0drms_a, r0drms_c;
    Eigen::Matrix<double, 3, 3> rot1_n, rot1_a, rot1_c, rot1_o;
    std::array<Eigen::Matrix<double, 3, MSAR>, 3> rot1_s;
    Eigen::Matrix<double, 3, 3> rot2_n, rot2_a, rot2_c, rot2_o;
    std::array<Eigen::Matrix<double, 3, MSAR>, 3> rot2_s;
    Eigen::Matrix<double, 3, 3> rot3_a, rot3_c, rot3_o;
    std::array<Eigen::Matrix<double, 3, MSAR>, 3> rot3_s;
    Eigen::Matrix<double, 3, 3> rot4_c;
    std::array<Eigen::Matrix<double, 3, MSAR>, 3> rot4_s;
    Eigen::Matrix<double, 3, 5> r1_n, r1_a, r1_c, r1_o;
    
    Eigen::Matrix<double, 3, 3> r_n_before, r_ca_before, r_c_before;
    Eigen::Matrix<double, 3, 3> r_n_after, r_ca_after, r_c_after;
    
    *n_soln = 0;
    t_ang[0] = pi;  // The omega dihedral angle, always has value 180 degrees
    t_ang[1] = pi;

    // First, we solve the loop closure problem for the current set of dihedrals, pre rotation
    initialize_loop_closure(b_len, b_ang, t_ang);
    solve_3pep_poly(r_n.col(1), r_a.col(1), r_a.col(3), r_c.col(3), r_soln_n, r_soln_a, r_soln_c, &n_soln_before);
    ctx->soln_no_before = n_soln_before;   
    if (ctx->soln_no_before == 0)
        ctx->soln_no_before = 1;
    r_n_before  = r_n.middleCols<3>(1);
    r_ca_before = r_a.middleCols<3>(1);
    r_c_before  = r_c.middleCols<3>(1);
    loop_Jacobian(r_n_before, r_ca_before, r_c_before, &before_jacobi);
    ctx->jacobi_before = before_jacobi;
    
    r0drms_n = r_n.middleCols<3>(1);
    r0drms_a = r_a.middleCols<3>(1);
    r0drms_c = r_c.middleCols<3>(1);

    // driver angle
    // We do a driver rotation
    // Note that the second argument of the function driver_rot is coordiate of
    // the target atom AFTER rotation, and this one is computed in the driver_rot function
    // See tripep_closure.h for what all these arguments mean

    if (ctx->mc.is_phi) { // We rotate the phi dihedral, which goes from N to CA
        driver_rot(r_c.col(3), r1_c.col(3), r_a.col(4), r_n.col(4), dih_ch);
        driver_rot(r_a.col(3), r1_a.col(3), r_a.col(4), r_n.col(4), dih_ch);
        r_c.col(3) = r1_c.col(3);
        r_a.col(3) = r1_a.col(3);
    } else {
        driver_rot(r_o.col(0), r1_o.col(0), r_a.col(0), r_c.col(0), dih_ch);
        driver_rot(r_n.col(1), r1_n.col(1), r_a.col(0), r_c.col(0), dih_ch);
        driver_rot(r_a.col(1), r1_a.col(1), r_a.col(0), r_c.col(0), dih_ch);
        r_o.col(0) = r1_o.col(0);
        r_n.col(1) = r1_n.col(1);
        r_a.col(1) = r1_a.col(1);
    }

    r0_n = r_n.middleCols<3>(1);
    r0_a = r_a.middleCols<3>(1);
    r0_c = r_c.middleCols<3>(1);
    r0_o = r_o.middleCols<3>(1);
    auto r0_s = std::span(r_s).subspan<1, 3>();

    solve_3pep_poly(r_n.col(1), r_a.col(1), r_a.col(3), r_c.col(3), r_soln_n, r_soln_a, r_soln_c, n_soln);  // Find new loop closure solution post rotation

    if (calc_rmsd) { // KP NOTE: Can be remvoed
        // for (k = 0; k < *n_soln; k++) {
        //     r_n.col(1) = r_soln_n[k].col(0);
        //     r_a.col(1) = r_soln_a[k].col(0);
        //     r_c.col(1) = r_soln_c[k].col(0);
        //     sum = 0.0e0;
        //     for (i = 0; i < 3; i++) {
        //         sum += (r_soln_n[k] - r0drms_n).squaredNorm();
        //         sum += (r_soln_a[k] - r0drms_a).squaredNorm();
        //         sum += (r_soln_c[k] - r0drms_c).squaredNorm();
        //     }
        //     rmsd = sqrt(sum / 9.0e0);
        //     if (rmsd < min_value) {
        //         min_value = rmsd;
        //         index     = k;
        //     }
        // }
        if (*n_soln > 0) {
            index = (int)(threefryrand() * (*n_soln));
            if (index >= (*n_soln)) {
                fprintf(sim->STATUS,
                        "ERROR!!!: index %d is greater than equal to no. of solutions %d.\n", index,
                        *n_soln);
                exit(1);
            }
        }
        if ((*n_soln > 0) && (*n_soln <= deg_pol)) {
            if (ctx->mc.is_phi) { // This is the case where driver angle is res. 4. (0~4)
                // rotation with n0-a0
                for (z = 0; z < ns[1]; z++)
                    get_rot(r0_s[0].col(z), r_soln_s[0].col(z), r_c.col(0), r_soln_n[index].col(0),
                            r_soln_a[index].col(0), r0_c.col(0), r_soln_c[index].col(0));
                get_rot(r0_o.col(0), rot1_o.col(0), r_c.col(0), r_soln_n[index].col(0), r_soln_a[index].col(0), r0_c.col(0),
                        r_soln_c[index].col(0));
                get_rot(r0_n.col(1), rot1_n.col(1), r_c.col(0), r_soln_n[index].col(0), r_soln_a[index].col(0), r0_c.col(0),
                        r_soln_c[index].col(0));
                get_rot(r0_a.col(1), rot1_a.col(1), r_c.col(0), r_soln_n[index].col(0), r_soln_a[index].col(0), r0_c.col(0),
                        r_soln_c[index].col(0));
                
                for (z = 0; z < ns[2]; z++)
                    get_rot(r0_s[1].col(z), rot1_s[1].col(z), r_c.col(0), r_soln_n[index].col(0),
                            r_soln_a[index].col(0), r0_c.col(0), r_soln_c[index].col(0));
                get_rot(r0_c.col(1), rot1_c.col(1), r_c.col(0), r_soln_n[index].col(0), r_soln_a[index].col(0), r0_c.col(0),
                        r_soln_c[index].col(0));

                // rotation with a0-c0
                get_rot(rot1_o.col(0), r_soln_o[index].col(0), r_soln_n[index].col(0), r_soln_a[index].col(0),
                        r_soln_c[index].col(0), rot1_n.col(1), r_soln_n[index].col(1));
                get_rot(rot1_a.col(1), rot2_a.col(1), r_soln_n[index].col(0), r_soln_a[index].col(0),
                        r_soln_c[index].col(0), rot1_n.col(1), r_soln_n[index].col(1));
                for (z = 0; z < ns[2]; z++)
                    get_rot(rot1_s[1].col(z), rot2_s[1].col(z), r_soln_n[index].col(0), r_soln_a[index].col(0),
                            r_soln_c[index].col(0), rot1_n.col(1), r_soln_n[index].col(1));
                get_rot(rot1_c.col(1), rot2_c.col(1), r_soln_n[index].col(0), r_soln_a[index].col(0),
                        r_soln_c[index].col(0), rot1_n.col(1), r_soln_n[index].col(1));

                // rotation with c0-n1
                for (z = 0; z < ns[2]; z++)
                    get_rot(rot2_s[1].col(z), rot3_s[1].col(z), r_soln_a[index].col(0), r_soln_c[index].col(0),
                            r_soln_n[index].col(1), rot2_a.col(1), r_soln_a[index].col(1));
                get_rot(rot2_c.col(1), rot3_c.col(1), r_soln_a[index].col(0), r_soln_c[index].col(0),
                        r_soln_n[index].col(1), rot2_a.col(1), r_soln_a[index].col(1));

                // rotation with n1-a1
                for (z = 0; z < ns[2]; z++)
                    get_rot(rot3_s[1].col(z), r_soln_s[1].col(z), r_soln_c[index].col(0), r_soln_n[index].col(1),
                            r_soln_a[index].col(1), rot3_c.col(1), r_soln_c[index].col(1));

                // rotation with a[+]-n[+]
                driver_rot(r0_o.col(2), r_soln_o[index].col(2), r_a.col(4), r_n.col(4), dih_ch);
                for (z = 0; z < ns[3]; z++)
                    driver_rot(r0_s[2].col(z), rot1_s[2].col(z), r_a.col(4), r_n.col(4), dih_ch);
                driver_rot(r0_n.col(2), rot1_n.col(2), r_a.col(4), r_n.col(4), dih_ch);
                driver_rot(r0_c.col(1), rot1_c.col(1), r_a.col(4), r_n.col(4), dih_ch);
                driver_rot(r0_o.col(1), rot1_o.col(1), r_a.col(4), r_n.col(4), dih_ch);
                driver_rot(r0_a.col(1), rot1_a.col(1), r_a.col(4), r_n.col(4), dih_ch);

                // rotation with c2-a2
                for (z = 0; z < ns[3]; z++)
                    get_rot(rot1_s[2].col(z), r_soln_s[2].col(z), r_n.col(4), r_soln_c[index].col(2),
                            r_soln_a[index].col(2), rot1_n.col(2), r_soln_n[index].col(2));
                get_rot(rot1_c.col(1), rot2_c.col(1), r_n.col(4), r_soln_c[index].col(2), r_soln_a[index].col(2),
                        rot1_n.col(2), r_soln_n[index].col(2));
                get_rot(rot1_o.col(1), rot2_o.col(1), r_n.col(4), r_soln_c[index].col(2), r_soln_a[index].col(2),
                        rot1_n.col(2), r_soln_n[index].col(2));
                get_rot(rot1_a.col(1), rot2_a.col(1), r_n.col(4), r_soln_c[index].col(2), r_soln_a[index].col(2),
                        rot1_n.col(2), r_soln_n[index].col(2));

                // rotation with a2-n2
                get_rot(rot2_o.col(1), rot3_o.col(1), r_soln_c[index].col(2), r_soln_a[index].col(2),
                        r_soln_n[index].col(2), rot2_c.col(1), r_soln_c[index].col(1));
                get_rot(rot2_a.col(1), rot3_a.col(1), r_soln_c[index].col(2), r_soln_a[index].col(2),
                        r_soln_n[index].col(2), rot2_c.col(1), r_soln_c[index].col(1));

                // rotation with n2-c1
                get_rot(rot3_o.col(1), r_soln_o[index].col(1), r_soln_a[index].col(2), r_soln_n[index].col(2),
                        r_soln_c[index].col(1), rot3_a.col(1), r_soln_a[index].col(1));
            } else
            // This is the case where driver angle is res. 0. (0~4)
            {
                // rotation with n[-]-a[-]
                for (z = 0; z < ns[1]; z++)
                    driver_rot(r0_s[0].col(z), rot1_s[0].col(z), r_a.col(0), r_c.col(0), dih_ch);
                driver_rot(r0_c.col(0), rot1_c.col(0), r_a.col(0), r_c.col(0), dih_ch);
                driver_rot(r0_o.col(0), rot1_o.col(0), r_a.col(0), r_c.col(0), dih_ch);
                driver_rot(r0_n.col(1), rot1_n.col(1), r_a.col(0), r_c.col(0), dih_ch);
                driver_rot(r0_a.col(1), rot1_a.col(1), r_a.col(0), r_c.col(0), dih_ch);
                for (z = 0; z < ns[2]; z++)
                    driver_rot(r0_s[1].col(z), rot1_s[1].col(z), r_a.col(0), r_c.col(0), dih_ch);
                driver_rot(r0_c.col(1), rot1_c.col(1), r_a.col(0), r_c.col(0), dih_ch);

                // rotation with n0-a0
                for (z = 0; z < ns[1]; z++)
                    get_rot(rot1_s[0].col(z), r_soln_s[0].col(z), r_c.col(0), r_soln_n[index].col(0),
                            r_soln_a[index].col(0), rot1_c.col(0), r_soln_c[index].col(0));
                get_rot(rot1_o.col(0), rot2_o.col(0), r_c.col(0), r_soln_n[index].col(0), r_soln_a[index].col(0),
                        rot1_c.col(0), r_soln_c[index].col(0));
                get_rot(rot1_n.col(1), rot2_n.col(1), r_c.col(0), r_soln_n[index].col(0), r_soln_a[index].col(0),
                        rot1_c.col(0), r_soln_c[index].col(0));
                get_rot(rot1_a.col(1), rot2_a.col(1), r_c.col(0), r_soln_n[index].col(0), r_soln_a[index].col(0),
                        rot1_c.col(0), r_soln_c[index].col(0));
                for (z = 0; z < ns[2]; z++)
                    get_rot(rot1_s[1].col(z), rot2_s[1].col(z), r_c.col(0), r_soln_n[index].col(0),
                            r_soln_a[index].col(0), rot1_c.col(0), r_soln_c[index].col(0));
                get_rot(rot1_c.col(1), rot2_c.col(1), r_c.col(0), r_soln_n[index].col(0), r_soln_a[index].col(0),
                        rot1_c.col(0), r_soln_c[index].col(0));

                // rotation with a0-c0
                get_rot(rot2_o.col(0), r_soln_o[index].col(0), r_soln_n[index].col(0), r_soln_a[index].col(0),
                        r_soln_c[index].col(0), rot2_n.col(1), r_soln_n[index].col(1));
                get_rot(rot2_a.col(1), rot3_a.col(1), r_soln_n[index].col(0), r_soln_a[index].col(0),
                        r_soln_c[index].col(0), rot2_n.col(1), r_soln_n[index].col(1));
                for (z = 0; z < ns[2]; z++)
                    get_rot(rot2_s[1].col(z), rot3_s[1].col(z), r_soln_n[index].col(0), r_soln_a[index].col(0),
                            r_soln_c[index].col(0), rot2_n.col(1), r_soln_n[index].col(1));
                get_rot(rot2_c.col(1), rot3_c.col(1), r_soln_n[index].col(0), r_soln_a[index].col(0),
                        r_soln_c[index].col(0), rot2_n.col(1), r_soln_n[index].col(1));

                // rotation with c0-n1
                for (z = 0; z < ns[2]; z++)
                    get_rot(rot3_s[1].col(z), rot4_s[1].col(z), r_soln_a[index].col(0), r_soln_c[index].col(0),
                            r_soln_n[index].col(1), rot3_a.col(1), r_soln_a[index].col(1));
                get_rot(rot3_c.col(1), rot4_c.col(1), r_soln_a[index].col(0), r_soln_c[index].col(0),
                        r_soln_n[index].col(1), rot3_a.col(1), r_soln_a[index].col(1));

                // rotation with n1-a1
                for (z = 0; z < ns[2]; z++)
                    get_rot(rot4_s[1].col(z), r_soln_s[1].col(z), r_soln_c[index].col(0), r_soln_n[index].col(1),
                            r_soln_a[index].col(1), rot4_c.col(1), r_soln_c[index].col(1));

                for (i = 0; i < 3; i++)
                    r_soln_o[index].col(2)[i] = r0_o.col(2)[i];

                // rotation with c2-a2
                for (z = 0; z < ns[3]; z++)
                    get_rot(r0_s[2].col(z), r_soln_s[2].col(z), r_n.col(4), r_soln_c[index].col(2),
                            r_soln_a[index].col(2), r0_n.col(2), r_soln_n[index].col(2));
                get_rot(r0_c.col(1), rot1_c.col(1), r_n.col(4), r_soln_c[index].col(2), r_soln_a[index].col(2), r0_n.col(2),
                        r_soln_n[index].col(2));
                get_rot(r0_o.col(1), rot1_o.col(1), r_n.col(4), r_soln_c[index].col(2), r_soln_a[index].col(2), r0_n.col(2),
                        r_soln_n[index].col(2));
                get_rot(r0_a.col(1), rot1_a.col(1), r_n.col(4), r_soln_c[index].col(2), r_soln_a[index].col(2), r0_n.col(2),
                        r_soln_n[index].col(2));

                // rotation with a2-n2
                get_rot(rot1_o.col(1), rot2_o.col(1), r_soln_c[index].col(2), r_soln_a[index].col(2),
                        r_soln_n[index].col(2), rot1_c.col(1), r_soln_c[index].col(1));
                get_rot(rot1_a.col(1), rot2_a.col(1), r_soln_c[index].col(2), r_soln_a[index].col(2),
                        r_soln_n[index].col(2), rot1_c.col(1), r_soln_c[index].col(1));

                // rotation with n2-c1
                get_rot(rot2_o.col(1), r_soln_o[index].col(1), r_soln_a[index].col(2), r_soln_n[index].col(2),
                        r_soln_c[index].col(1), rot2_a.col(1), r_soln_a[index].col(1));
            }
            
            r_n_after = r_soln_n[index];
            r_ca_after = r_soln_a[index];
            r_c_after = r_soln_c[index];

            r_n.middleCols<3>(1) = r_soln_n[index];
            r_a.middleCols<3>(1) = r_soln_a[index];
            r_c.middleCols<3>(1) = r_soln_c[index];
            r_o.middleCols<3>(1) = r_soln_o[index];
            std::ranges::copy(std::span(r_s).subspan<1, 3>(), r_soln_s.begin() + 1);
        }
    }

    loop_Jacobian(r_n_after, r_ca_after, r_c_after, &after_jacobi);
    ctx->jacobi_after = after_jacobi;

    return;
}

//========================================================================================================

void get_template(
    struct Topology *top,
    struct Simulation *sim
) {
    FILE *the_file;

    char line[250], temp[10];
    int  initial_template[MAX_RES], mod_template[MAX_RES], len_template[MAX_RES];
    int  i, j;
    int  prev_tmpl = 0;
    int  len_tmpl  = 0;

    for (i = 0; i < top->nresidues; i++) {
        initial_template[i] = 0;
        len_template[i]     = 0;
    }

    if ((the_file = fopen(sim->template_file.c_str(), "r")) == NULL) {
        fprintf(sim->STATUS, "ERROR: Can't open the file: %s!\n", sim->template_file);
        exit(1);
    }
    while (fgets(line, 100, the_file) != NULL) {
        if (strncmp(line, "ATOM", 4) == 0) {
            strncpy(temp, &(line[22]), 4);
            strcpy(&temp[4], "\0");
            initial_template[atoi(temp) - 1] = 1;
        }
    }
    fclose(the_file);

    for (i = 0; i < top->nresidues; i++) {
        if ((prev_tmpl == 0) && (initial_template[i] == 1))
            len_tmpl = 1;
        else if ((prev_tmpl == 1) && (initial_template[i] == 1))
            len_tmpl++;
        else if ((prev_tmpl == 1) && (initial_template[i] == 0)) {
            for (j = i - len_tmpl; j < i; j++)
                len_template[j] = len_tmpl;
            len_tmpl = 0;
        }
        if ((i == top->nresidues - 1) && (initial_template[i] == 1))
            for (j = i - len_tmpl + 1; j <= i; j++)
                len_template[j] = len_tmpl;
        prev_tmpl = initial_template[i];
    }

    for (i = 0; i < top->nresidues; i++)
        if (len_template[i] < MIN_LEN_TEMPLATE)
            mod_template[i] = 0;
        else
            mod_template[i] = initial_template[i];

    prev_tmpl = 0;
    for (i = 0; i < top->nresidues; i++) {
        if ((prev_tmpl == 0) && (mod_template[i] == 1) && (i != 0))
            top->is_template[i] = 0;
        else if ((prev_tmpl == 1) && (mod_template[i] == 0))
            top->is_template[i - 1] = 0;
        else
            top->is_template[i] = mod_template[i];
        prev_tmpl = mod_template[i];
    }
    for (i = 0; i < top->nresidues; i++) {
        fprintf(sim->STATUS, "%5d%5d%5d%5d%5d\n", i, initial_template[i], len_template[i],
                mod_template[i], top->is_template[i]);
        fflush(sim->STATUS);
    }
    return;
}

//========================================================================================================
void check_bb(
    struct Topology *top, struct Context *ctx
) {
    Vec3 tmp1, tmp2, tmp3, tmp4, plane1, plane2, bisect1, bisect2;
    int           i;

    for (i = 0; i < top->nresidues - 2; i++) {
        MakeVector(ctx->native[ctx->native_residue[i].CA].xyz, ctx->native[ctx->native_residue[i].N].xyz, &tmp1);
        MakeVector(ctx->native[ctx->native_residue[i].CA].xyz, ctx->native[ctx->native_residue[i].O].xyz, &tmp2);
        MakeVector(ctx->native[ctx->native_residue[i + 2].CA].xyz, ctx->native[ctx->native_residue[i + 2].N].xyz,
                   &tmp3);
        MakeVector(ctx->native[ctx->native_residue[i + 2].CA].xyz, ctx->native[ctx->native_residue[i + 2].O].xyz,
                   &tmp4);
        CrossProduct(tmp1, tmp2, &plane1);
        CrossProduct(tmp3, tmp4, &plane2);
        ctx->a_PCA[i] = Angle(plane1, plane2);
        ctx->a_PCA[i] *= rad2deg;
        Normalize(&tmp1);
        Normalize(&tmp2);
        Normalize(&tmp3);
        Normalize(&tmp4);
        bisect(tmp1, tmp2, &bisect1);
        bisect(tmp3, tmp4, &bisect2);
        ctx->a_bCA[i] = Angle(bisect1, bisect2);
        ctx->a_bCA[i] *= rad2deg;

        MakeVector(ctx->native[ctx->native_residue[i].C].xyz, ctx->native[ctx->native_residue[i + 1].N].xyz, &tmp1);
        MakeVector(ctx->native[ctx->native_residue[i + 1].N].xyz, ctx->native[ctx->native_residue[i + 1].CA].xyz,
                   &tmp2);
        MakeVector(ctx->native[ctx->native_residue[i + 1].CA].xyz, ctx->native[ctx->native_residue[i + 1].C].xyz,
                   &tmp3);
        ctx->phim[i] = struct_calc_dih_ang(tmp1, tmp2, tmp3);
        ctx->phim[i] *= rad2deg;
        MakeVector(ctx->native[ctx->native_residue[i + 1].N].xyz, ctx->native[ctx->native_residue[i + 1].CA].xyz,
                   &tmp1);
        MakeVector(ctx->native[ctx->native_residue[i + 1].CA].xyz, ctx->native[ctx->native_residue[i + 1].C].xyz,
                   &tmp2);
        MakeVector(ctx->native[ctx->native_residue[i + 1].C].xyz, ctx->native[ctx->native_residue[i + 2].N].xyz, &tmp3);
        ctx->psim[i] = struct_calc_dih_ang(tmp1, tmp2, tmp3);
        ctx->psim[i] *= rad2deg;
        ctx->phim[i] += 180.;
        ctx->psim[i] += 180.;

        //    fprintf(STATUS, "%4d %4d %9.3f %9.3f %9.3f %9.3f %9.3f\n", i+1, i+2,  dih_CA[i],
        //    ang_CA[i], ang_CA[i+1], len12_CA[i], len03_CA[i]);
    }

    return;
}

//========================================================================================================
void initialize_torsion(
    struct Topology *top, struct System *sys, struct Simulation *sim
) {
    int r, cur_res, value;
    int x, y, z, w;

    FILE *ftor;
    char  line[1000];

    for (r = 0; r < top->nresidues; r++)
        for (x = 0; x < 6; x++)
            for (y = 0; y < 6; y++)
                for (z = 0; z < 6; z++)
                    for (w = 0; w < 6; w++)
                        sys->torsion_E[r][x][y][z][w] = 1000;

    fprintf(sim->STATUS, "Opening the file: %s\n", sim->triplet_file);
    if ((ftor = fopen(sim->triplet_file.c_str(), "r")) == NULL) {
        fprintf(sim->STATUS, "ERROR: Can't open the file: %s!\n", sim->triplet_file);
        exit(1);
    }
    while (fgets(line, 200, ftor) != NULL) {
        sscanf(line, "%d %*d%*d%*d %d %d %d %d %d %*s%*s", &cur_res, &x, &y, &z, &w, &value);
        sys->torsion_E[cur_res][x][y][z][w] = value;
    }
    fclose(ftor);

    return;
}

//========================================================================================================
void initialize_sct(
    struct Topology *top, struct System *sys, struct Simulation *sim
) {
    int r, cur_res, value;
    int x, y, z, w;
    printf("DEBUG: initialize_sct\n");
    FILE *ftor;
    char  line[1000];

    for (r = 0; r < top->nresidues; r++)
        for (x = 0; x < 12; x++)
            for (y = 0; y < 12; y++)
                for (z = 0; z < 12; z++)
                    for (w = 0; w < 12; w++)
                        sys->sct_E[r][x][y][z][w] = 1000;

    fprintf(sim->STATUS, "Opening the file: %s\n", sim->sctorsion_file);
    if ((ftor = fopen(sim->sctorsion_file.c_str(), "r")) == NULL) {
        fprintf(sim->STATUS, "ERROR: Can't open the file: %s!\n", sim->sctorsion_file);
        exit(1);
    }
    while (fgets(line, 200, ftor) != NULL) {
        sscanf(line, "%d %*d%*d%*d %d %d %d %d %d %*s%*s", &cur_res, &x, &y, &z, &w, &value);
        sys->sct_E[cur_res][x][y][z][w] = value;
        //    fprintf(sim->STATUS, "%10d\n", sys->sct_E[cur_res][x][y][z][w]);
    }
    fclose(ftor);
    printf("DEBUG: initialize_sct done\n");
    return;
}

//========================================================================================================
void initialize_aromatic(
    struct Topology *top, struct System *sys, struct Simulation *sim,
    struct Context *ctx, struct MCIntegrator *integrator
) {
    int value;
    int x, r;

    FILE *ftor;
    char  line[1000];

    fprintf(sim->STATUS, "Opening the file: %s\n", sim->aromatic_file.c_str());
    if ((ftor = fopen(sim->aromatic_file.c_str(), "r")) == NULL) {
        fprintf(sim->STATUS, "ERROR: Can't open the file: %s!\n", sim->aromatic_file.c_str());
        exit(1);
    }
    while (fgets(line, 200, ftor) != NULL) {
        sscanf(line, "%d %d %*s%*s", &x, &value);
        sys->aromatic_E[x] = value;
        //    fprintf(sim->STATUS, "%2d %5d\n", x, sys->aromatic_E[x]);
    }
    fclose(ftor);

    for (r = 0; r < top->nresidues; r++)
        if ((strcmp(ctx->native_residue[r].res, "PHE") == 0) ||
            (strcmp(ctx->native_residue[r].res, "TRP") == 0) ||
            (strcmp(ctx->native_residue[r].res, "RING") == 0)) {
            top->Res_aromatic[top->Naromatic] = r;
            //      fprintf(sim->STATUS, "%5d %5d\n", Naromatic, Res_aromatic[Naromatic]);
            top->Naromatic++;
        }
    return;
}



//========================================================================================================
void read_cluster(
    std::string filename, struct System *sys
) {
    std::ifstream file(filename);
    std::string line;
    
    // Skip header line if it exists
    std::getline(file, line);

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string value;
        
        int res, clu;
        float phi, psi;

        std::getline(ss, value, ','); clu = std::stoi(value);
        std::getline(ss, value, ','); res = std::stoi(value);
        std::getline(ss, value, ','); phi = std::stof(value);
        std::getline(ss, value, ','); psi = std::stof(value);
    
        sys->cluster_phi[res][clu] = phi;
        sys->cluster_psi[res][clu] = psi;
    }
    std::cout << "Cluster data read from " << filename << std::endl;
    file.close();
    //  for (i=0;i<20;i++)
    //    for (j=0;j<NOCLUSTERS;j++)
    //      fprintf(sim->STATUS, "%d %d %8.5f %8.5f\n", i, j, sys->cluster_phi[i][j], sys->cluster_psi[i][j]);
    return;
}

//========================================================================================================
void check_phipsi(
    const struct Topology *top, struct Context *ctx
) {
    Eigen::Matrix<double, 3, MAXSEQUENCE> res_CA;
    Eigen::Matrix<double, 3, MAXSEQUENCE> res_N;
    Eigen::Matrix<double, 3, MAXSEQUENCE> res_C;

    for (int i = 0; i < top->nresidues; i++) {
        res_N.col(i)  = ctx->native[ctx->native_residue[i].N].xyz;
        res_CA.col(i) = ctx->native[ctx->native_residue[i].CA].xyz;
        res_C.col(i)  = ctx->native[ctx->native_residue[i].C].xyz;
    }

    for (int i = 0; i < top->nresidues - 2; i++) {
        c_dih_ang(res_C.col(i), res_N.col(i + 1), res_CA.col(i + 1), res_C.col(i + 1), &ctx->cur_phi[i]);
        c_dih_ang(res_N.col(i + 1), res_CA.col(i + 1), res_C.col(i + 1), res_N.col(i + 2), &ctx->cur_psi[i]);
        ctx->cur_phi[i] *= rad2deg;
        ctx->cur_psi[i] *= rad2deg;
        //    fprintf(sim->STATUS, "%4d %4d %9.3f %9.3f %9.3f %9.3f %9.3f\n", i+1, i+2,  dih_CA[i],
        //    ang_CA[i], ang_CA[i+1], len12_CA[i], len03_CA[i]);
    }

    return;
}

//========================================================================================================
void initialize_secstr(
    struct Simulation *sim, struct Topology *top
) {
    FILE *ftor;
    char  line[1000];
    char  iscorrect[600], secstr_org[600];
    char  tmp_str[10];
    int   i;

    fprintf(sim->STATUS, "Opening the file: %s\n", sim->sec_str_file);
    if ((ftor = fopen(sim->sec_str_file.c_str(), "r")) == NULL) {
        fprintf(sim->STATUS, "ERROR: Can't open the file: %s!\n", sim->sec_str_file);
        exit(1);
    }
    fgets(line, 500, ftor);
    sscanf(line, "%s", iscorrect);
    fgets(line, 500, ftor);
    sscanf(line, "%s", secstr_org);
    fclose(ftor);

    for (i = 0; i < top->nresidues; i++) {
        sprintf(tmp_str, "%c", iscorrect[i]);
        if ((atoi(tmp_str) >= CUT_SECSTR) && secstr_org[i] == 'H')
            top->secstr[i] = secstr_org[i];
        else if ((atoi(tmp_str) >= CUT_SECSTR) && secstr_org[i] == 'E')
            top->secstr[i] = secstr_org[i];
        else if ((atoi(tmp_str) >= CUT_SECSTR) && secstr_org[i] == 'C')
            top->secstr[i] = 'L';
        else
            top->secstr[i] = 'C';
    }
    fprintf(sim->STATUS, "\n");
    // fprintf(sim->STATUS, "secondary structure correct?\n%s\n", iscorrect);
    // fprintf(sim->STATUS, "input secondary structure:\n%s\n", secstr_org);
    fprintf(sim->STATUS, "secondary structure:\n%s\n", top->secstr);

    return;
}