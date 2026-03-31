#ifndef LOOP_H
#define LOOP_H

#include "globals.h"
#include "in_out.h"

void integloop(
    struct Context *ctx, struct Topology *top, struct MCIntegrator *integrator, struct System *sys, struct Simulation *sim,
    float step_size, int *n_soln);
void get_orgco(
    struct Topology *top, struct Context *ctx,
    double rorg_n[5][3], double rorg_a[5][3], double rorg_c[5][3], double rorg_o[5][3],
    double rorg_s[5][MSAR][3], int ns[5], int n0);
void get_orgat(
    struct Context *ctx,
    double rorg_n[5][3], double rorg_a[5][3], double rorg_c[5][3], double rorg_o[5][3],
    double rorg_s[5][MSAR][3], int ns[5], int r, int a);
void get_coord(
    struct Topology *top, struct Context *ctx,
    double r_n[5][3], double r_a[5][3], double r_c[5][3], double r_o[5][3],
    double r_s[5][MSAR][3], int ns[5], int n0);
void get_atom(
    struct Context *ctx,
    double r_n[5][3], double r_a[5][3], double r_c[5][3], double r_o[5][3],
    double r_s[5][MSAR][3], int ns[5], int r, int i);
void put_coord(
    struct Topology *top, struct Context *ctx,
    double r_n[5][3], double r_a[5][3], double r_c[5][3], double r_o[5][3],
    double r_s[5][MSAR][3], int ns[5], int n0);
void put_atom(struct Context *ctx,
              double r_n[5][3], double r_a[5][3], double r_c[5][3], double r_o[5][3],
              double r_s[5][MSAR][3], int ns[5], int r, int a);
void yangloop(
    struct Context *ctx, struct Simulation *sim,
    double r_n[5][3], double r_a[5][3], double r_c[5][3], double r_o[5][3],
    double r_s[5][MSAR][3], int ns[5], double b_len[6], double b_ang[7], int n0,
    float dih_ch, char res_name[5][4], int *n_soln);
void yang_rotate(
    struct Topology *top, struct MCIntegrator *integrator, struct Context *ctx,
    int n0, int is_phi);
void get_template(
    struct Topology *top,
    struct Simulation *sim
);
void check_bb(
    struct Topology *top, struct Context *ctx
);
void check_phipsi(
    const struct Topology *top, struct Context *ctx
);
void initialize_torsion(
    struct Topology *top, struct System *sys, struct Simulation *sim
);
void initialize_sct(
    struct Topology *top, struct System *sys, struct Simulation *sim
);
void initialize_aromatic(
    struct Topology *top, struct System *sys, struct Simulation *sim,
    struct Context *ctx, struct MCIntegrator *integrator
);
int  will2edo(int);
void read_cluster(
    struct Simulation *sim, struct System *sys
);
void initialize_secstr(
    struct Simulation *sim, struct Topology *top
);
void Read_movable_region(
    struct Simulation *sim,
    struct Topology *top
);

#endif