#ifndef TRIPEP_CLOSURE_H
#define TRIPEP_CLOSURE_H

#include "globals.h"
using namespace MCPU2;

#define max_soln 16
#define deg_pol 16
#include "sturm.h"

#define max(a, b) ((a) > (b)) ? (a) : (b)
#define min(a, b) ((a) < (b)) ? (a) : (b)

Vec3::Scalar dot_product(const Vec3& va, const Vec3& vb);
void   matmul(const Mat3& ma, const Vec3& mb, Vec3& mc);
double sign(double a, double b);
void   solve_3pep_poly(const Vec3& r_n1, const Vec3& r_a1, const Vec3& r_a3, const Vec3& r_c3,
                     std::array<Eigen::Matrix3d, max_soln>& r_soln_n,
                     std::array<Eigen::Matrix3d, max_soln>& r_soln_a,
                     std::array<Eigen::Matrix3d, max_soln>& r_soln_c, int *n_soln) ;
void   initialize_loop_closure(double b_len[6], double b_ang[7], double t_ang[2]);
void get_input_angles(int *n_soln, const Vec3& r_n1, const Vec3& r_a1, const Vec3& r_a3, const Vec3& r_c3);
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
void coord_from_poly_roots(int *n_soln, double roots[max_soln], const Vec3& r_n1, const Vec3& r_a1,
                           const Vec3& r_a3, const Vec3& r_c3,
                           std::array<Eigen::Matrix3d, max_soln>& r_soln_n,
                           std::array<Eigen::Matrix3d, max_soln> &r_soln_a,
                           std::array<Eigen::Matrix3d, max_soln>& r_soln_c);
double calc_t2(double t0);
double calc_t1(double t0, double t2);

void   calc_dih_ang(const Vec3& r1, const Vec3& r2, const Vec3& r3, double *angle);
double  struct_calc_dih_ang(const Vec3& r1, const Vec3& r2, const Vec3& r3);
void   calc_bnd_ang(const Vec3& r1, const Vec3& r2, double *angle);
void   c_bnd_len(const Vec3& r1, const Vec3& r2, double *length);
void   c_bnd_ang(const Vec3& r1, const Vec3& r2, const Vec3& r3, double *angle);
void   c_dih_ang(const Vec3& r1, const Vec3& r2, const Vec3& r3, const Vec3& r4, double *angle);
void cross(const Vec3& p, const Vec3& q, Vec3& s);
void quaternion(const Vec3& axis, double quarter_ang, Eigen::Quaterniond& p);
void rotation_matrix(const Eigen::Quaterniond& q, Eigen::Matrix3d& U);
void apply_rotation_eigen(const Eigen::Ref<Vec3> obj_i, 
                          Eigen::Ref<Vec3> obj_f, 
                          const Eigen::Ref<Vec3> tor2, 
                          const Eigen::Ref<Vec3> tor3, 
                          double total_angle);
void get_rot(const Eigen::Ref<Vec3> obj_i, Eigen::Ref<Vec3> obj_f, 
             const Eigen::Ref<Vec3> tor1, const Eigen::Ref<Vec3> tor2, const Eigen::Ref<Vec3> tor3,
             const Eigen::Ref<Vec3> tor4_i, const Eigen::Ref<Vec3> tor4_f);
void driver_rot(const Eigen::Ref<Vec3> obj_i, Eigen::Ref<Vec3> obj_f, 
                const Eigen::Ref<Vec3> tor2, const Eigen::Ref<Vec3> tor3,
                double dih_change);

#endif