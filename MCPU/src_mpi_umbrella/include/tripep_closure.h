#ifndef TRIPEP_CLOSURE_H
#define TRIPEP_CLOSURE_H

#define max_soln 16
#define deg_pol 16
#include "sturm.h"

#define max(a, b) ((a) > (b)) ? (a) : (b)
#define min(a, b) ((a) < (b)) ? (a) : (b)

double dot_product(double va[3], double vb[3]);
void   matmul(double ma[3][3], double mb[3], double mc[3]);
double sign(double a, double b);
void   solve_3pep_poly(double r_n1[3], double r_a1[3], double r_a3[3], double r_c3[3],
                       double r_soln_n[max_soln][3][3], double r_soln_a[max_soln][3][3],
                       double r_soln_c[max_soln][3][3], int *n_soln);
void   initialize_loop_closure(double b_len[6], double b_ang[7], double t_ang[2]);
void get_input_angles(int *n_soln, double r_n1[3], double r_a1[3], double r_a3[3], double r_c3[3]);
void test_two_cone_existence_soln(double tt, double kx, double et, double ap, int *n_soln,
                                  char cone_type[2]);
void get_poly_coeff(double poly_coeff[deg_pol + 1]);
void poly_mul_sub2(double u1[5][5], double u2[5][5], double u3[5][5], double u4[5][5], int p1[2],
                   int p2[2], int p3[2], int p4[2], double u5[5][5], int p5[2]);
void poly_mul2(double u1[5][5], double u2[5][5], int p1[2], int p2[2], double u3[5][5], int p3[2]);
void poly_sub2(double u1[5][5], double u2[5][5], int p1[2], int p2[2], double u3[5][5], int p3[2]);
void poly_mul_sub1(double u1[17], double u2[17], double u3[17], double u4[17], int p1, int p2,
                   int p3, int p4, double u5[17], int *p5);
void poly_mul1(double u1[17], double u2[17], int p1, int p2, double u3[17], int *p3);
void poly_sub1(double u1[17], double u2[17], int p1, int p2, double u3[17], int *p3);
void coord_from_poly_roots(int *n_soln, double roots[max_soln], double r_n1[3], double r_a1[3],
                           double r_a3[3], double r_c3[3], double r_soln_n[max_soln][3][3],
                           double r_soln_a[max_soln][3][3], double r_soln_c[max_soln][3][3]);
double calc_t2(double t0);
double calc_t1(double t0, double t2);
void   calc_dih_ang(double r1[3], double r2[3], double r3[3], double *angle);
void   calc_bnd_ang(double r1[3], double r2[3], double *angle);
void   c_bnd_len(double r1[3], double r2[3], double *length);
void   c_bnd_ang(double r1[3], double r2[3], double r3[3], double *angle);
void   c_dih_ang(double r1[3], double r2[3], double r3[3], double r4[3], double *angle);
float  struct_calc_dih_ang(struct vector r1, struct vector r2, struct vector r3);
void   cross(double p[3], double q[3], double s[3]);
void   quaternion(double axis[3], double quarter_ang, double p[4]);
void   rotation_matrix(double q[4], double U[3][3]);
void   get_rot(double obj_i[3], double obj_f[3], double tor1[3], double tor2[3], double tor3[3],
               double tor4_i[3], double tor_f[3]);
void driver_rot(double obj_i[3], double obj_f[3], double tor2[3], double tor3[3], float dih_change);

#endif