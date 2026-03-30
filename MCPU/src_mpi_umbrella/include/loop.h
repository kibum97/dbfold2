#ifndef LOOP_H
#define LOOP_H

#include "in_out.h"

void integloop(float step_size, int *n_soln);
void get_orgco(double rorg_n[5][3], double rorg_a[5][3], double rorg_c[5][3], double rorg_o[5][3],
               double rorg_s[5][MSAR][3], int ns[5], int n0);
void get_orgat(double rorg_n[5][3], double rorg_a[5][3], double rorg_c[5][3], double rorg_o[5][3],
               double rorg_s[5][MSAR][3], int ns[5], int r, int i);
void get_coord(double r_n[5][3], double r_a[5][3], double r_c[5][3], double r_o[5][3],
               double r_s[5][MSAR][3], int ns[5], int n0);
void get_atom(double r_n[5][3], double r_a[5][3], double r_c[5][3], double r_o[5][3],
              double r_s[5][MSAR][3], int ns[5], int r, int i);
void put_coord(double r_n[5][3], double r_a[5][3], double r_c[5][3], double r_o[5][3],
               double r_s[5][MSAR][3], int ns[5], int n0);
void put_atom(double r_n[5][3], double r_a[5][3], double r_c[5][3], double r_o[5][3],
              double r_s[5][MSAR][3], int ns[5], int r, int i);
void yangloop(double r_n[5][3], double r_a[5][3], double r_c[5][3], double r_o[5][3],
              double r_s[5][MSAR][3], int ns[5], double b_len[6], double b_ang[7], int n0,
              float dih_ch, char res_name[5][4], int *n_soln);
void yang_rotate(int n0, int is_phi);
void get_template();
void check_bb();
void check_phipsi();
void initialize_torsion();
void initialize_sct();
void initialize_aromatic();
int  will2edo(int);
void read_cluster();
void initialize_secstr();
void Read_movable_region();

#endif