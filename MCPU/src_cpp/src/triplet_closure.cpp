#include "tripep_closure.h"
#include "vector.h"

//!----------------------------------------------------------------------
//! Copyright (C) 2003
//!      Chaok Seok, Evangelos Coutsias, Matthew Jacobson, and Ken Dill
//!      UCSF and Univeristy of New Mexico
//! Witten by Chaok Seok 2003.
//!----------------------------------------------------------------------------
//!----------------------------------------------------------------------------
//!*************** Tripeptide Loop Closure Algorithm *****************
//!----------------------------------------------------------------------------
// MODULE tripep_closure
//!----------------------------------------------------------------------------
#define pi 3.141592653589793238462643383279502884197e0
#define two_pi 2.0e0 * pi
#define deg2rad pi / 180.0e0
#define rad2deg 180.0e0 / pi
#define max(a, b) ((a) > (b)) ? (a) : (b)
#define min(a, b) ((a) < (b)) ? (a) : (b)
int print_level = 0;
double len0[6], b_ang0[7], t_ang0[2];
double aa13_min_sqr, aa13_max_sqr;
double delta[4], xi[3], eta[3], alpha[3], theta[3];
double cos_alpha[3], sin_alpha[3], cos_theta[3], sin_theta[3];
double cos_delta[4], sin_delta[4];
double cos_xi[3], cos_eta[3], sin_xi[3], sin_eta[3];
Vec3 r_a1a3, r_a1n1, r_a3c3;
double b_a1a3[3], b_a1n1[3], b_a3c3[3];
double len_na[3], len_ac[3], len_aa[3];
double C0[3][3], C1[3][3], C2[3][3];
double Q[5][17], R[3][17];

Vec3::Scalar dot_product(const Vec3& va, const Vec3& vb) {
    return va.dot(vb);
}

void matmul(const Mat3& ma, const Vec3& mb, Vec3& mc) {
    mc = ma * mb;
}

double sign(double a, double b) {
    if (b >= 0.)
        return fabs(a);
    else
        return -fabs(a);
}

//!-----------------------------------------------------------------------
void solve_3pep_poly(const Vec3& r_n1, const Vec3& r_a1, const Vec3& r_a3, const Vec3& r_c3,
                     std::array<Eigen::Matrix3d, max_soln>& r_soln_n,
                     std::array<Eigen::Matrix3d, max_soln>& r_soln_a,
                     std::array<Eigen::Matrix3d, max_soln>& r_soln_c, int *n_soln) {

    double poly_coeff[deg_pol + 1], roots[max_soln];
    int    var_deg_pol;

    //  call get_input_angles(n_soln, r_n1, r_a1, r_a3, r_c3)
    get_input_angles(n_soln, r_n1, r_a1, r_a3, r_c3);

    //  if (n_soln == 0) then
    //     return
    //  end if
    if (*n_soln == 0)
        return;

    //  call get_poly_coeff(poly_coeff)
    get_poly_coeff(poly_coeff);

    //  call solve_sturm(deg_pol, n_soln, poly_coeff, roots)
    var_deg_pol = deg_pol;
    solve_sturm(&var_deg_pol, n_soln, poly_coeff, roots);

    //  if (n_soln == 0) then
    //!     print*, 'return 2'
    //     return
    //  end if
    if (*n_soln == 0)
        return;

    //  call coord_from_poly_roots(n_soln, roots, r_n1, r_a1, r_a3, r_c3, r_soln_n, r_soln_a,
    //  r_soln_c)
    coord_from_poly_roots(n_soln, roots, r_n1, r_a1, r_a3, r_c3, r_soln_n, r_soln_a, r_soln_c);

    return;
    // end subroutine solv_3pep_poly
}
//!-----------------------------------------------------------------------
// subroutine initialize_loop_closure(b_len, b_ang, t_ang)
void initialize_loop_closure(double b_len[6], double b_ang[7], double t_ang[2]) {
    //!-----------------------------------------------------------------------
    //! Input angles for the given bond lengths and angles
    //! i.e. the fixed values! Namely the fixed omega dihedral angle and known bond lengths and
    //! atomic angles based on basic chemistry, orbitals etc
    //!-----------------------------------------------------------------------
    //  implicit none
    double len1, len2, a_min, a_max;
    Vec3 rr_a1, rr_c1, rr_n2, rr_a2, rr_n2a2_ref, rr_c1a1;
    Vec3 rr_a1a2, dr, bb_c1a1, bb_a1a2, bb_a2n2;
    Eigen::Quaterniond p;
    Eigen::Matrix3d Us;
    Vec3 mulpro;
    Vec3 tmp_val;
    //  real(dp), parameter :: tol_secant = 1.0d-15
    double tol_secant = 1.0e-15;
    //  integer, parameter :: max_iter_sturm = 100, max_iter_secant = 20
    int max_iter_sturm  = 100;
    int max_iter_secant = 20;
    //  integer :: i
    int i, j;

    //  call initialize_sturm(tol_secant, max_iter_sturm, max_iter_secant)
    initialize_sturm(&tol_secant, &max_iter_sturm, &max_iter_secant);

    //  len0(1:6) = b_len(1:6)
    std::copy(b_len, b_len + 6, len0);
    std::copy(b_ang, b_ang + 7, b_ang0);
    std::copy(t_ang, t_ang + 2, t_ang0);

    //  rr_c1(1:3) = 0.0d0
    rr_c1.setZero();
    Vec3 axis{1.0, 0.0, 0.0};

    //  do i = 0, 1
    for (i = 0; i < 2; i++)  // Loop through all atoms in triplet of residues
    {
        double s1, c1, s2, c2;
        sincos(b_ang0[3 * i + 1], &s1, &c1);
        sincos(b_ang0[3 * i + 2], &s2, &c2);

        double l0 = len0[3 * i];
        double l1 = len0[3 * i + 1];
        double l2 = len0[3 * i + 2];

        rr_a1 = Vec3{c1 * l0, s1 * l0, 0.0};
        rr_n2 = Vec3{l1, 0.0e0, 0.0e0};  // LEFT OFF HERE!
        rr_c1a1 = rr_a1 - rr_c1;
        rr_n2a2_ref = Vec3{-c2 * l2, s2 * l2, 0.0e0};
        
        
        quaternion(axis, t_ang0[i] * 0.25e0, p);
        rotation_matrix(p, Us);
        matmul(Us, rr_n2a2_ref, mulpro);
        
        rr_a2 = mulpro + rr_n2;
        rr_a1a2 = rr_a2 - rr_a1;
        dr = rr_a1a2;

        len2 = dot_product(dr, dr);
        len1 = sqrt(len2);
        len_aa[i + 1] = len1;

        bb_c1a1 = rr_c1a1 / l0;
        bb_a1a2 = rr_a1a2 / len1;
        bb_a2n2 = (rr_n2 - rr_a2) / l2;

        tmp_val -= bb_a1a2;

        calc_bnd_ang(tmp_val, bb_a2n2, &xi[i + 1]);
        //     ! eta
        //     call calc_bnd_ang(bb_a1a2, -bb_c1a1, eta(i+1))
        for (j = 0; j < 3; j++)
            tmp_val[j] = -bb_c1a1[j];
        calc_bnd_ang(bb_a1a2, tmp_val, &eta[i]);
        //     ! delta: pi -  dih of N(1)CA(1)CA(3)C(3)
        //     call calc_dih_ang(bb_c1a1, bb_a1a2, bb_a2n2, delta(i+1))
        calc_dih_ang(bb_c1a1, bb_a1a2, bb_a2n2, &delta[i + 1]);
        //     delta(i+1) = pi - delta(i+1)
        delta[i + 1] = pi - delta[i + 1];
    }
    //  end do

    //  a_min = b_ang(4) - (xi(2) + eta(2))
    a_min = b_ang[3] - (xi[1] + eta[1]);
    //  a_max = min(b_ang(4) + (xi(2) + eta(2)), pi)
    a_max = min((b_ang[3] + (xi[1] + eta[1])), pi);

    //  ! min/max of base length
    //!  print*, 'len1, len3=', len_aa(2:3)
    //  fprintf(STATUS,"len1, len3= %9.5f %9.5f\n", len_aa[1], len_aa[2]);
    //!  print*, 'a_min, a_max=', a_min*rad2deg, a_max*rad2deg
    //  fprintf(STATUS,"a_min, a_max= %9.5f %9.5f\n", a_min*rad2deg, a_max*rad2deg);
    //  aa13_min_sqr = len_aa(2)**2 + len_aa(3)**2 - 2.0d0*len_aa(2)*len_aa(3)*cos(a_min)
    aa13_min_sqr =
        pow(len_aa[1], 2) + pow(len_aa[2], 2) - 2.0e0 * len_aa[1] * len_aa[2] * cos(a_min);
    //  aa13_max_sqr = len_aa(2)**2 + len_aa(3)**2 - 2.0d0*len_aa(2)*len_aa(3)*cos(a_max)
    aa13_max_sqr =
        pow(len_aa[1], 2) + pow(len_aa[2], 2) - 2.0e0 * len_aa[1] * len_aa[2] * cos(a_max);
    //!  print*, 'aa13_min_sqr,aa13_max_sqr', aa13_min_sqr,aa13_max_sqr
    //  fprintf(STATUS,"aa13_min_sqr,aa13_max_sqr %9.5f %9.5f\n", aa13_min_sqr, aa13_max_sqr);

    // end subroutine initialize_loop_closure
}
//!-----------------------------------------------------------------------
// subroutine get_input_angles(n_soln, r_n1, r_a1, r_a3, r_c3)
void get_input_angles(int *n_soln, const Vec3& r_n1, const Vec3& r_a1, const Vec3& r_a3, const Vec3& r_c3) {
    //!-----------------------------------------------------------------------
    //! Input angles and vectors (later used in coordinates)
    //!-----------------------------------------------------------------------
    Vec3 tmp_val;
    int i;
    char cone_type[2];

    *n_soln = max_soln;

    Vec3 r_a1a3 = r_a3 - r_a1;
    double dr_sqr = r_a1a3.squaredNorm();
    len_aa[0] = sqrt(dr_sqr);

    if ((dr_sqr < aa13_min_sqr) || (dr_sqr > aa13_max_sqr)) {
        *n_soln = 0;
        return;
    }

    Vec3 r_a1n1 = r_n1 - r_a1;
    len_na[0] = r_a1n1.norm();
    len_na[1] = len0[2];
    len_na[2] = len0[5];
    Vec3 r_a3c3 = r_c3 - r_a3;
    len_ac[0] = len0[0];
    len_ac[1] = len0[3];
    len_ac[2] = r_a3c3.norm();

    //  ! unit vectors
    Vec3 b_a1n1 = r_a1n1 / len_na[0];
    Vec3 b_a3c3 = r_a3c3 / len_ac[2];
    Vec3 b_a1a3 = r_a1a3 / len_aa[0];


    tmp_val = -b_a1n1;
    calc_dih_ang(tmp_val, b_a1a3, b_a3c3, &delta[3]);
    delta[0] = delta[3];

    tmp_val = -b_a1a3;
    calc_bnd_ang(tmp_val, b_a1n1, &xi[0]);

    calc_bnd_ang(b_a1a3, b_a3c3, &eta[2]);

    // KP - STOPED HERE - NEED TO FINISH THIS FUNCTION
    //  do i = 1, 3
    for (i = 0; i < 3; i++) {
        //     cos_delta(i) = cos(delta(i))
        cos_delta[i + 1] = cos(delta[i + 1]);
        //     sin_delta(i) = sin(delta(i))
        sin_delta[i + 1] = sin(delta[i + 1]);
        //     cos_xi(i) = cos(xi(i))
        //     cos_xi(i) = cos(xi(i))
        cos_xi[i] = cos(xi[i]);
        //     sin_xi(i) = sin(xi(i))
        sin_xi[i] = sin(xi[i]);
        //     sin_xi(i) = sin(xi(i))
        sin_xi[i] = sin(xi[i]);
        //     cos_eta(i) = cos(eta(i))
        cos_eta[i] = cos(eta[i]);
        //     cos_eta(i) = cos(eta(i))
        cos_eta[i] = cos(eta[i]);
        //     sin_eta(i) = sin(eta(i))
        sin_eta[i] = sin(eta[i]);
        //     sin_eta(i) = sin(eta(i))
        sin_eta[i] = sin(eta[i]);
        //  end do
    }
    //  cos_delta(0) = cos_delta(3)
    cos_delta[0] = cos_delta[3];
    //  sin_delta(0) = sin_delta(3)
    sin_delta[0] = sin_delta[3];

    //  ! theta (N, CA, C) bond angle
    //  theta(1) = b_ang0(1)
    theta[0] = b_ang0[0];
    //  theta(2) = b_ang0(4)
    theta[1] = b_ang0[3];
    //  theta(3) = b_ang0(7)
    theta[2] = b_ang0[6];
    //  do i = 1, 3
    //     cos_theta(i) = cos(theta(i))
    //  end do
    for (i = 0; i < 3; i++)
        cos_theta[i] = cos(theta[i]);

    //  ! alpha
    //  cos_alpha(1) = -(len_aa(1)**2 + len_aa(2)**2 - len_aa(3)**2)/(2.0d0*len_aa(1)*len_aa(2))
    cos_alpha[0] = -(pow(len_aa[0], 2) + pow(len_aa[1], 2) - pow(len_aa[2], 2)) /
                   (2.0e0 * len_aa[0] * len_aa[1]);
    //  alpha(1) = acos(cos_alpha(1))
    alpha[0] = acos(cos_alpha[0]);
    //  sin_alpha(1) = sin(alpha(1))
    sin_alpha[0] = sin(alpha[0]);
    //  cos_alpha(2) = (len_aa(2)**2 + len_aa(3)**2 - len_aa(1)**2)/(2.0d0*len_aa(2)*len_aa(3))
    cos_alpha[1] = (pow(len_aa[1], 2) + pow(len_aa[2], 2) - pow(len_aa[0], 2)) /
                   (2.0e0 * len_aa[1] * len_aa[2]);
    //  alpha(2) = acos(cos_alpha(2))
    alpha[1] = acos(cos_alpha[1]);
    //  sin_alpha(2) = sin(alpha(2))
    sin_alpha[1] = sin(alpha[1]);
    //  alpha(3) = pi - alpha(1) + alpha(2)
    alpha[2] = pi - alpha[0] + alpha[1];
    //  cos_alpha(3) = cos(alpha(3))
    cos_alpha[2] = cos(alpha[2]);
    //  sin_alpha(3) = sin(alpha(3))
    sin_alpha[2] = sin(alpha[2]);

    //  if (print_level > 0) then
    // if (print_level > 0) {
    //     //     write(*,'(a,3f9.4)') 'xi = ', xi(1:3)*rad2deg
    //     fprintf(STATUS, "xi = %9.4f%9.4f%9.4f\n", xi[0] * rad2deg, xi[1] * rad2deg,
    //             xi[2] * rad2deg);
    //     //     write(*,'(a,3f9.4)') 'eta = ', eta(1:3)*rad2deg
    //     fprintf(STATUS, "eta = %9.4f%9.4f%9.4f\n", eta[0] * rad2deg, eta[1] * rad2deg,
    //             eta[2] * rad2deg);
    //     //     write(*,'(a,3f9.4)') 'delta = ', delta(1:3)*rad2deg
    //     fprintf(STATUS, "delta = %9.4f%9.4f%9.4f\n", delta[1] * rad2deg, delta[2] * rad2deg,
    //             delta[3] * rad2deg);
    //     //     write(*,'(a,3f9.4)') 'theta = ', theta(1:3)*rad2deg
    //     fprintf(STATUS, "theta = %9.4f%9.4f%9.4f\n", theta[0] * rad2deg, theta[1] * rad2deg,
    //             theta[2] * rad2deg);
    //     //     write(*,'(a,3f9.4)') 'alpha = ', alpha(1:3)*rad2deg
    //     fprintf(STATUS, "alpha = %9.4f%9.4f%9.4f\n", alpha[0] * rad2deg, alpha[1] * rad2deg,
    //             alpha[2] * rad2deg);
    //     //  end if
    // }

    //  ! check for existence of soln
    //  do i = 1, 3
    for (i = 0; i < 3; i++) {
        //     call test_two_cone_existence_soln(theta(i), xi(i), eta(i), alpha(i), &
        //          n_soln, cone_type)
        test_two_cone_existence_soln(theta[i], xi[i], eta[i], alpha[i], n_soln, cone_type);
        //     if (n_soln == 0) then
        //        print*, 'return 1', i
        //        return
        //     end if
        if (*n_soln == 0)
            return;
        //  end do
    }

    return;
    // end subroutine get_input_angles
}
//!-----------------------------------------------------------------------
// subroutine test_two_cone_existence_soln(tt, kx, et, ap, n_soln, cone_type)
void test_two_cone_existence_soln(double tt, double kx, double et, double ap, int *n_soln,
                                  char cone_type[2]) {
    //  implicit none
    //  real(dp), intent(in) :: tt, kx, et, ap
    //  integer, intent(out) :: n_soln
    //  character(len=2), intent(out) :: cone_type
    //  character(len=2) :: case_type
    //  real(dp) :: at, ex, abs_at, ap1, kx1, et1
    double at, ex, abs_at, ap1, kx1, et1;
    //  real(dp) :: cos_tx1, cos_tx2, cos_te1, cos_te2, cos_ea1, cos_ea2, cos_xa1, cos_xa2
    double cos_tx1, cos_tx2, cos_te1, cos_te2, cos_ea1, cos_ea2, cos_xa1, cos_xa2;
    //  logical :: s1, s2, t1, t2, complicated = .false.
    int s1, s2, t1, t2;
    int complicated = 0;
    //  real(dp), parameter :: half_pi = 0.5d0*pi

    //  n_soln = max_soln
    *n_soln = max_soln;

    //  ap1 = ap
    ap1 = ap;
    //  kx1 = kx
    kx1 = kx;
    //  et1 = et
    et1 = et;

    //  at = ap1 - tt
    at = ap1 - tt;
    //  ex = kx1 + et1
    ex = kx1 + et1;
    //  abs_at = abs(at)
    abs_at = fabs(at);

    //  ! case of no soln
    //  if (abs_at > ex) then
    //     n_soln = 0
    //     return
    //  end if
    if (abs_at > ex) {
        *n_soln = 0;
        return;
    }

    //  if (complicated) then
    //     ! find type of intersection
    //     cos_tx1 = cos(tt+kx1)
    //     cos_tx2 = cos(tt-kx1)
    //     cos_te1 = cos(tt+et1)
    //     cos_te2 = cos(tt-et1)
    //     cos_ea1 = cos(et1+ap1)
    //     cos_ea2 = cos(et1-ap1)
    //     cos_xa1 = cos(kx1+ap1)
    //     cos_xa2 = cos(kx1-ap1)
    //     s1 = .false.; s2 = .false.; t1 = .false.; t2 = .false.
    //     if ((cos_te1-cos_xa2)*(cos_te1-cos_xa1) <= 0.0d0) s1 = .true.
    //     if ((cos_te2-cos_xa2)*(cos_te2-cos_xa1) <= 0.0d0) s2 = .true.
    //     if ((cos_tx1-cos_ea2)*(cos_tx1-cos_ea1) <= 0.0d0) t1 = .true.
    //     if ((cos_tx2-cos_ea2)*(cos_tx2-cos_ea1) <= 0.0d0) t2 = .true.
    //  end if
    if (complicated) {
        //    cos_tx1 = cos(tt+kx1)
        cos_tx1 = cos(tt + kx1);
        //    cos_tx2 = cos(tt-kx1)
        cos_tx2 = cos(tt - kx1);
        //    cos_te1 = cos(tt+et1)
        cos_te1 = cos(tt + et1);
        //    cos_te2 = cos(tt-et1)
        cos_te2 = cos(tt - et1);
        //    cos_ea1 = cos(et1+ap1)
        cos_ea1 = cos(et1 + ap1);
        //    cos_ea2 = cos(et1-ap1)
        cos_ea2 = cos(et1 - ap1);
        //    cos_xa1 = cos(kx1+ap1)
        cos_xa1 = cos(kx1 + ap1);
        //    cos_xa2 = cos(kx1-ap1)
        cos_xa2 = cos(kx1 - ap1);
        //    s1 = .false.; s2 = .false.; t1 = .false.; t2 = .false.
        s1 = 0;
        s2 = 0;
        t1 = 0;
        t2 = 0;
        //    if ((cos_te1-cos_xa2)*(cos_te1-cos_xa1) <= 0.0d0) s1 = .true.
        if ((cos_te1 - cos_xa2) * (cos_te1 - cos_xa1) <= 0.0e0)
            s1 = 0;
        //    if ((cos_te2-cos_xa2)*(cos_te2-cos_xa1) <= 0.0d0) s2 = .true.
        if ((cos_te2 - cos_xa2) * (cos_te2 - cos_xa1) <= 0.0e0)
            s2 = 0;
        //    if ((cos_tx1-cos_ea2)*(cos_tx1-cos_ea1) <= 0.0d0) t1 = .true.
        if ((cos_tx1 - cos_ea2) * (cos_tx1 - cos_ea1) <= 0.0e0)
            t1 = 0;
        //    if ((cos_tx2-cos_ea2)*(cos_tx2-cos_ea1) <= 0.0d0) t2 = .true.
        if ((cos_tx2 - cos_ea2) * (cos_tx2 - cos_ea1) <= 0.0e0)
            t2 = 0;
    }

    return;
    // end subroutine test_two_cone_existence_soln
}
//!-----------------------------------------------------------------------
// subroutine get_poly_coeff(poly_coeff)
void get_poly_coeff(double poly_coeff[deg_pol + 1]) {
    //  implicit none
    //  real(dp), intent(out) :: poly_coeff(0:deg_pol)
    //  integer :: i, j
    int i, j;
    //  real(dp) :: A0, A1, A2, A3, A4, A21, A22, A31, A32, A41, A42
    double A0, A1, A2, A3, A4, A21, A22, A31, A32, A41, A42;
    //  real(dp) :: B0(3), B1(3), B2(3), B3(3), B4(3), B5(3), B6(3), B7(3), B8(3)
    double B0[3], B1[3], B2[3], B3[3], B4[3], B5[3], B6[3], B7[3], B8[3];
    //  real(dp), dimension(0:4,0:4) :: u11, u12, u13, u31, u32, u33
    double u11[5][5], u12[5][5], u13[5][5], u31[5][5], u32[5][5], u33[5][5];
    //  real(dp), dimension(0:4,0:4) :: um1, um2, um3, um4, um5, um6, q_tmp
    double um1[5][5], um2[5][5], um3[5][5], um4[5][5], um5[5][5], um6[5][5], q_tmp[5][5];
    //  integer, dimension(2) :: p1, p3, p_um1, p_um2, p_um3, p_um4, p_um5, p_um6, p_Q
    int p1[2], p3[2], p_um1[2], p_um2[2], p_um3[2], p_um4[2], p_um5[2], p_um6[2], p_Q[2];
    //  integer :: p2, p4, p_f1, p_f2, p_f3, p_f4, p_f5, p_f6, p_f7, &
    //       p_f8, p_f9, p_f10, p_f11, p_f12, p_f13, p_f14, p_f15, p_f16, p_f17, &
    //       p_f18, p_f19, p_f20, p_f21, p_f22, p_f23, p_f24, p_f25, p_f26, p_f27
    int p2, p4, p_f1, p_f2, p_f3, p_f4, p_f5, p_f6, p_f7, p_f8, p_f9;
    int p_f10, p_f11, p_f12, p_f13, p_f14, p_f15, p_f16, p_f17, p_f18;
    int p_f19, p_f20, p_f21, p_f22, p_f23, p_f24, p_f25, p_f26;
    //  integer :: p_final
    int p_final;
    //  real(dp), dimension(0:16) :: f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, &
    //       f12, f13, f14, f15, f16, f17, f18, f19, f20, f21, f22, f23, f24, f25, f26, f27
    double f1[17], f2[17], f3[17], f4[17], f5[17], f6[17], f7[17], f8[17], f9[17];
    double f10[17], f11[17], f12[17], f13[17], f14[17], f15[17], f16[17], f17[17], f18[17];
    double f19[17], f20[17], f21[17], f22[17], f23[17], f24[17], f25[17], f26[17];

    //  ! A0, B0
    //  do i = 1, 3
    for (i = 0; i < 3; i++) {
        //     A0 = cos_alpha(i)*cos_xi(i)*cos_eta(i) - cos_theta(i)
        A0 = cos_alpha[i] * cos_xi[i] * cos_eta[i] - cos_theta[i];
        //     A1 = -sin_alpha(i)*cos_xi(i)*sin_eta(i)
        A1 = -sin_alpha[i] * cos_xi[i] * sin_eta[i];
        //     A2 = sin_alpha(i)*sin_xi(i)*cos_eta(i)
        A2 = sin_alpha[i] * sin_xi[i] * cos_eta[i];
        //     A3 = sin_xi(i)*sin_eta(i)
        A3 = sin_xi[i] * sin_eta[i];
        //     A4 = A3*cos_alpha(i)
        A4 = A3 * cos_alpha[i];
        //     j = i - 1
        j = i;
        //     A21 = A2*cos_delta(j)
        A21 = A2 * cos_delta[j];
        //     A22 = A2*sin_delta(j)
        A22 = A2 * sin_delta[j];
        //     A31 = A3*cos_delta(j)
        A31 = A3 * cos_delta[j];
        //     A32 = A3*sin_delta(j)
        A32 = A3 * sin_delta[j];
        //     A41 = A4*cos_delta(j)
        A41 = A4 * cos_delta[j];
        //     A42 = A4*sin_delta(j)
        A42 = A4 * sin_delta[j];
        //     B0(i) = A0 + A22 + A31
        B0[i] = A0 + A22 + A31;
        //     B1(i) = 2.0d0*(A1 + A42)
        B1[i] = 2.0e0 * (A1 + A42);
        //     B2(i) = 2.0d0*(A32 - A21)
        B2[i] = 2.0e0 * (A32 - A21);
        //     B3(i) = -4.0d0*A41
        B3[i] = -4.0e0 * A41;
        //     B4(i) = A0 + A22 - A31
        B4[i] = A0 + A22 - A31;
        //     B5(i) = A0 - A22 - A31
        B5[i] = A0 - A22 - A31;
        //     B6(i) = -2.0d0*(A21 + A32)
        B6[i] = -2.0e0 * (A21 + A32);
        //     B7(i) = 2.0d0*(A1 - A42)
        B7[i] = 2.0e0 * (A1 - A42);
        //     B8(i) = A0 - A22 + A31
        B8[i] = A0 - A22 + A31;
        //  end do
    }

    //  ! C0i
    //  i = 1
    i = 0;
    //  C0(0:2,i) = (/ B0(i), B2(i), B5(i) /);
    C0[i][0] = B0[i];
    C0[i][1] = B2[i];
    C0[i][2] = B5[i];
    //  C1(0:2,i) = (/ B1(i), B3(i), B7(i) /);
    C1[i][0] = B1[i];
    C1[i][1] = B3[i];
    C1[i][2] = B7[i];
    //  C2(0:2,i) = (/ B4(i), B6(i), B8(i) /);
    C2[i][0] = B4[i];
    C2[i][1] = B6[i];
    C2[i][2] = B8[i];
    //  do i = 2, 3
    for (i = 1; i < 3; i++) {
        //     C0(0:2,i) = (/ B0(i), B1(i), B4(i) /)
        C0[i][0] = B0[i];
        C0[i][1] = B1[i];
        C0[i][2] = B4[i];
        //     C1(0:2,i) = (/ B2(i), B3(i), B6(i) /)
        C1[i][0] = B2[i];
        C1[i][1] = B3[i];
        C1[i][2] = B6[i];
        //     C2(0:2,i) = (/ B5(i), B7(i), B8(i) /)
        C2[i][0] = B5[i];
        C2[i][1] = B7[i];
        C2[i][2] = B8[i];
        //  end do
    }

    //  ! first determinant
    //  do i = 0, 2
    for (i = 0; i < 3; i++) {
        //     u11(i,0) = C0(i,1)
        u11[0][i] = C0[0][i];
        //     u12(i,0) = C1(i,1)
        u12[0][i] = C1[0][i];
        //     u13(i,0) = C2(i,1)
        u13[0][i] = C2[0][i];
        //     u31(0,i) = C0(i,2) STRANGE !!!
        u31[i][0] = C0[1][i];
        //     u32(0,i) = C1(i,2) STRANGE !!!
        u32[i][0] = C1[1][i];
        //     u33(0,i) = C2(i,2) STRANGE !!!
        u33[i][0] = C2[1][i];
        //  end do
    }

    //  p1(1:2) = (/ 2, 0 /)
    p1[0] = 2;
    p1[1] = 0;
    //  p3(1:2) = (/ 0, 2 /)
    p3[0] = 0;
    p3[1] = 2;

    //  call poly_mul_sub2(u32, u32, u31, u33, p3, p3, p3, p3, um1, p_um1)
    poly_mul_sub2(u32, u32, u31, u33, p3, p3, p3, p3, um1, p_um1);
    //  call poly_mul_sub2(u12, u32, u11, u33, p1, p3, p1, p3, um2, p_um2)
    poly_mul_sub2(u12, u32, u11, u33, p1, p3, p1, p3, um2, p_um2);
    //  call poly_mul_sub2(u12, u33, u13, u32, p1, p3, p1, p3, um3, p_um3)
    poly_mul_sub2(u12, u33, u13, u32, p1, p3, p1, p3, um3, p_um3);
    //  call poly_mul_sub2(u11, u33, u31, u13, p1, p3, p3, p1, um4, p_um4)
    poly_mul_sub2(u11, u33, u31, u13, p1, p3, p3, p1, um4, p_um4);
    //  call poly_mul_sub2(u13, um1, u33, um2, p1, p_um1, p3, p_um2, um5, p_um5)
    poly_mul_sub2(u13, um1, u33, um2, p1, p_um1, p3, p_um2, um5, p_um5);
    //  call poly_mul_sub2(u13, um4, u12, um3, p1, p_um4, p1, p_um3, um6, p_um6)
    poly_mul_sub2(u13, um4, u12, um3, p1, p_um4, p1, p_um3, um6, p_um6);
    //  call poly_mul_sub2(u11, um5, u31, um6, p1, p_um5, p3, p_um6, q_tmp, p_Q)
    poly_mul_sub2(u11, um5, u31, um6, p1, p_um5, p3, p_um6, q_tmp, p_Q);

    //  Q(0:4,0:4) = q_tmp(0:4,0:4)
    for (i = 0; i < 5; i++)
        for (j = 0; j < 5; j++)
            Q[i][j] = q_tmp[i][j];

    //  ! second determinant
    //  R(:,:) = 0.0d0
    for (i = 0; i < 3; i++)
        for (j = 0; j < 17; j++)
            R[i][j] = 0.;
    //  R(0:2,0) = C0(0:2,3)
    //  R(0:2,1) = C1(0:2,3)
    //  R(0:2,2) = C2(0:2,3)
    for (i = 0; i < 3; i++) {
        R[0][i] = C0[2][i];
        R[1][i] = C1[2][i];
        R[2][i] = C2[2][i];
    }
    //  p2 = 2
    p2 = 2;
    //  p4 = 4
    p4 = 4;

    //  call poly_mul_sub1(R(:,1), R(:,1), R(:,0), R(:,2), p2, p2, p2, p2, f1, p_f1)
    poly_mul_sub1(R[1], R[1], R[0], R[2], p2, p2, p2, p2, f1, &p_f1);
    //  call poly_mul1(R(:,1), R(:,2), p2, p2, f2, p_f2)
    poly_mul1(R[1], R[2], p2, p2, f2, &p_f2);
    //  call poly_mul_sub1(R(:,1), f1, R(:,0), f2, p2, p_f1, p2, p_f2, f3, p_f3)
    poly_mul_sub1(R[1], f1, R[0], f2, p2, p_f1, p2, p_f2, f3, &p_f3);
    //  call poly_mul1(R(:,2), f1, p2, p_f1, f4, p_f4)
    poly_mul1(R[2], f1, p2, p_f1, f4, &p_f4);
    //  call poly_mul_sub1(R(:,1), f3, R(:,0), f4, p2, p_f3, p2, p_f4, f5, p_f5)
    poly_mul_sub1(R[1], f3, R[0], f4, p2, p_f3, p2, p_f4, f5, &p_f5);

    //  call poly_mul_sub1(Q(:,1), R(:,1), Q(:,0), R(:,2), p4, p2, p4, p2, f6, p_f6)
    poly_mul_sub1(Q[1], R[1], Q[0], R[2], p4, p2, p4, p2, f6, &p_f6);
    //  call poly_mul_sub1(Q(:,2), f1, R(:,2), f6, p4, p_f1, p2, p_f6, f7, p_f7)
    poly_mul_sub1(Q[2], f1, R[2], f6, p4, p_f1, p2, p_f6, f7, &p_f7);
    //  call poly_mul_sub1(Q(:,3), f3, R(:,2), f7, p4, p_f3, p2, p_f7, f8, p_f8)
    poly_mul_sub1(Q[3], f3, R[2], f7, p4, p_f3, p2, p_f7, f8, &p_f8);
    //  call poly_mul_sub1(Q(:,4), f5, R(:,2), f8, p4, p_f5, p2, p_f8, f9, p_f9)
    poly_mul_sub1(Q[4], f5, R[2], f8, p4, p_f5, p2, p_f8, f9, &p_f9);

    //  call poly_mul_sub1(Q(:,3), R(:,1), Q(:,4), R(:,0), p4, p2, p4, p2, f10, p_f10)
    poly_mul_sub1(Q[3], R[1], Q[4], R[0], p4, p2, p4, p2, f10, &p_f10);
    //  call poly_mul_sub1(Q(:,2), f1, R(:,0), f10, p4, p_f1, p2, p_f10, f11, p_f11)
    poly_mul_sub1(Q[2], f1, R[0], f10, p4, p_f1, p2, p_f10, f11, &p_f11);
    //  call poly_mul_sub1(Q(:,1), f3, R(:,0), f11, p4, p_f3, p2, p_f11, f12, p_f12)
    poly_mul_sub1(Q[1], f3, R[0], f11, p4, p_f3, p2, p_f11, f12, &p_f12);

    //  call poly_mul_sub1(Q(:,2), R(:,1), Q(:,1), R(:,2), p4, p2, p4, p2, f13, p_f13)
    poly_mul_sub1(Q[2], R[1], Q[1], R[2], p4, p2, p4, p2, f13, &p_f13);
    //  call poly_mul_sub1(Q(:,3), f1, R(:,2), f13, p4, p_f1, p2, p_f13, f14, p_f14)
    poly_mul_sub1(Q[3], f1, R[2], f13, p4, p_f1, p2, p_f13, f14, &p_f14);
    //  call poly_mul_sub1(Q(:,3), R(:,1), Q(:,2), R(:,2), p4, p2, p4, p2, f15, p_f15)
    poly_mul_sub1(Q[3], R[1], Q[2], R[2], p4, p2, p4, p2, f15, &p_f15);
    //  call poly_mul_sub1(Q(:,4), f1, R(:,2), f15, p4, p_f1, p2, p_f15, f16, p_f16)
    poly_mul_sub1(Q[4], f1, R[2], f15, p4, p_f1, p2, p_f15, f16, &p_f16);
    //  call poly_mul_sub1(Q(:,1), f14, Q(:,0), f16, p4, p_f14, p4, p_f16, f17, p_f17)
    poly_mul_sub1(Q[1], f14, Q[0], f16, p4, p_f14, p4, p_f16, f17, &p_f17);

    //  call poly_mul_sub1(Q(:,2), R(:,2), Q(:,3), R(:,1), p4, p2, p4, p2, f18, p_f18)
    poly_mul_sub1(Q[2], R[2], Q[3], R[1], p4, p2, p4, p2, f18, &p_f18);
    //  call poly_mul_sub1(Q(:,1), R(:,2), Q(:,3), R(:,0), p4, p2, p4, p2, f19, p_f19)
    poly_mul_sub1(Q[1], R[2], Q[3], R[0], p4, p2, p4, p2, f19, &p_f19);
    //  call poly_mul_sub1(Q(:,3), f19, Q(:,2), f18, p4, p_f19, p4, p_f18, f20, p_f20)
    poly_mul_sub1(Q[3], f19, Q[2], f18, p4, p_f19, p4, p_f18, f20, &p_f20);
    //  call poly_mul_sub1(Q(:,1), R(:,1), Q(:,2), R(:,0), p4, p2, p4, p2, f21, p_f21)
    poly_mul_sub1(Q[1], R[1], Q[2], R[0], p4, p2, p4, p2, f21, &p_f21);
    //  call poly_mul1(Q(:,4), f21, p4, p_f21, f22, p_f22)
    poly_mul1(Q[4], f21, p4, p_f21, f22, &p_f22);
    //  call poly_sub1(f20, f22, p_f20, p_f22, f23, p_f23)
    poly_sub1(f20, f22, p_f20, p_f22, f23, &p_f23);
    //  call poly_mul1(R(:,0), f23, p2, p_f23, f24, p_f24)
    poly_mul1(R[0], f23, p2, p_f23, f24, &p_f24);
    //  call poly_sub1(f17, f24, p_f17, p_f24, f25, p_f25)
    poly_sub1(f17, f24, p_f17, p_f24, f25, &p_f25);
    //  call poly_mul_sub1(Q(:,4), f12, R(:,2), f25, p4, p_f12, p2, p_f25, f26, p_f26)
    poly_mul_sub1(Q[4], f12, R[2], f25, p4, p_f12, p2, p_f25, f26, &p_f26);
    //  call poly_mul_sub1(Q(:,0), f9, R(:,0), f26, p4, p_f9, p2, p_f26, poly_coeff, p_final)
    poly_mul_sub1(Q[0], f9, R[0], f26, p4, p_f9, p2, p_f26, poly_coeff, &p_final);

    //  if (p_final /= 16) then
    //     print*, 'Error. Degree of polynomial is not 16!'
    //     stop
    //  end if
    // if (p_final != 16) {
    //     fprintf(STATUS, "Error. Degree of polynomial is not 16!\n");
    //     exit(1);
    // }

    //  if (poly_coeff(16) < 0.0d0) then
    //     poly_coeff(0:16) = -poly_coeff(0:16)
    //  end if
    if (poly_coeff[16] < 0.0e0)
        for (i = 0; i < 17; i++)
            poly_coeff[i] *= -1.0;

    //  if (print_level > 0) then
    //     print*, 'poly_coeff'
    //     do i = 0, 16
    //        write(*,"(i5,e15.6)") i, poly_coeff(i)
    //     end do
    //  end if
    // if (print_level > 0) {
    //     fprintf(STATUS, "poly_coeff\n");
    //     for (i = 0; i < 17; i++)
    //         fprintf(STATUS, "%5d%15.6f\n", i, poly_coeff[i]);
    // }

    return;
    // end subroutine get_poly_coeff
}
//!----------------------------------------------------------------------------
// subroutine poly_mul_sub2(u1, u2, u3, u4, p1, p2, p3, p4, u5, p5)
void poly_mul_sub2(double u1[5][5], double u2[5][5], double u3[5][5], double u4[5][5], int p1[2],
                   int p2[2], int p3[2], int p4[2], double u5[5][5], int p5[2]) {
    //  implicit none
    //  real(dp), dimension(0:4,0:4), intent(in) :: u1, u2, u3, u4
    //  integer, dimension(2), intent(in) :: p1, p2, p3, p4
    //  real(dp), dimension(0:4,0:4), intent(out) :: u5
    //  integer, dimension(2), intent(out) :: p5
    //  real(dp), dimension(0:4,0:4) :: d1, d2
    double d1[5][5], d2[5][5];
    //  integer, dimension(2) :: pd1, pd2
    int pd1[2], pd2[2];

    //  call poly_mul2(u1, u2, p1, p2, d1, pd1)
    poly_mul2(u1, u2, p1, p2, d1, pd1);
    //  call poly_mul2(u3, u4, p3, p4, d2, pd2)
    poly_mul2(u3, u4, p3, p4, d2, pd2);
    //  call poly_sub2(d1, d2, pd1, pd2, u5, p5)
    poly_sub2(d1, d2, pd1, pd2, u5, p5);

    // end subroutine poly_mul_sub2
}
//!----------------------------------------------------------------------------
// subroutine poly_mul2(u1, u2, p1, p2, u3, p3)
void poly_mul2(double u1[5][5], double u2[5][5], int p1[2], int p2[2], double u3[5][5], int p3[2]) {
    //  implicit none
    //  real(dp), dimension(0:4,0:4), intent(in) :: u1, u2
    //  integer, dimension(2), intent(in) :: p1, p2
    //  real(dp), dimension(0:4,0:4), intent(out) :: u3
    //  integer, intent(out) :: p3(2)
    //  integer :: i1, j1, i2, j2, i3, j3, p11, p12, p21, p22
    int i1, j1, i2, j2, i3, j3, p11, p12, p21, p22;
    int i, j;
    //  real(dp) :: u1ij
    double u1ij;

    //  p3(:) = p1(:) + p2(:)
    for (i = 0; i < 2; i++)
        p3[i] = p1[i] + p2[i];
    for (i = 0; i < 5; i++)
        for (j = 0; j < 5; j++)
            u3[i][j] = 0.0e0;

    //  p11 = p1(1)
    p11 = p1[0];
    //  p12 = p1(2)
    p12 = p1[1];
    //  p21 = p2(1)
    p21 = p2[0];
    //  p22 = p2(2)
    p22 = p2[1];

    //  do i1 = 0, p12
    for (i1 = 0; i1 <= p12; i1++) {
        //     do j1 = 0, p11
        for (j1 = 0; j1 <= p11; j1++) {
            //        u1ij = u1(j1,i1)
            u1ij = u1[i1][j1];
            //        do i2 = 0, p22
            for (i2 = 0; i2 <= p22; i2++) {
                //           i3 = i1 + i2
                i3 = i1 + i2;
                //           do j2 = 0, p21
                for (j2 = 0; j2 <= p21; j2++) {
                    //              j3 = j1 + j2
                    j3 = j1 + j2;
                    //              u3(j3,i3) = u3(j3,i3) + u1ij*u2(j2,i2)
                    u3[i3][j3] = u3[i3][j3] + u1ij * u2[i2][j2];
                    //           end do
                }
                //        end do
            }
            //     end do
        }
        //  end do
    }

    // end subroutine poly_mul2
}
//!----------------------------------------------------------------------------
// subroutine poly_sub2(u1, u2, p1, p2, u3, p3)
void poly_sub2(double u1[5][5], double u2[5][5], int p1[2], int p2[2], double u3[5][5], int p3[2]) {
    //  implicit none
    //  real(dp), dimension(0:4,0:4), intent(in) :: u1, u2
    //  integer, intent(in) :: p1(2), p2(2)
    //  real(dp), dimension(0:4,0:4), intent(out) :: u3
    //  integer, intent(out) :: p3(2)
    //  integer :: i, j, p11, p12, p21, p22
    int i, j, p11, p12, p21, p22;
    //  logical :: i1_ok, i2_ok
    int i1_ok, i2_ok;

    //  p11 = p1(1)
    p11 = p1[0];
    //  p12 = p1(2)
    p12 = p1[1];
    //  p21 = p2(1)
    p21 = p2[0];
    //  p22 = p2(2)
    p22 = p2[1];
    //  p3(1) = max(p11,p21)
    p3[0] = max(p11, p21);
    //  p3(2) = max(p12,p22)
    p3[1] = max(p12, p22);

    //  do i = 0, p3(2)
    for (i = 0; i <= p3[1]; i++) {
        //     i1_ok = (i > p12)
        i1_ok = (i > p12);
        //     i2_ok = (i > p22)
        i2_ok = (i > p22);
        //     do j = 0, p3(1)
        for (j = 0; j <= p3[0]; j++) {
            //        if (i2_ok .or. (j > p21)) then
            //           u3(j,i) = u1(j,i)
            if (i2_ok || (j > p21))
                u3[i][j] = u1[i][j];
            //        else if (i1_ok .or. (j > p11)) then
            //           u3(j,i) = -u2(j,i)
            else if (i1_ok || (j > p11))
                u3[i][j] = -u2[i][j];
            //        else
            //           u3(j,i) = u1(j,i) - u2(j,i)
            else
                u3[i][j] = u1[i][j] - u2[i][j];
            //        end if
            //     end do
        }
        //  end do
    }

    return;
    // end subroutine poly_sub2
}
//!----------------------------------------------------------------------------
// subroutine poly_mul_sub1(u1, u2, u3, u4, p1, p2, p3, p4, u5, p5)
void poly_mul_sub1(double u1[17], double u2[17], double u3[17], double u4[17], int p1, int p2,
                   int p3, int p4, double u5[17], int *p5) {
    //  implicit none
    //  real(dp), dimension(0:16), intent(in) :: u1, u2, u3, u4
    //  integer, intent(in) :: p1, p2, p3, p4
    //  real(dp), dimension(0:16), intent(out) :: u5
    //  integer, intent(out) :: p5
    //  real(dp), dimension(0:16) :: d1, d2
    double d1[17], d2[17];
    //  integer :: pd1, pd2
    int pd1, pd2;

    //  call poly_mul1(u1, u2, p1, p2, d1, pd1)
    poly_mul1(u1, u2, p1, p2, d1, &pd1);
    //  call poly_mul1(u3, u4, p3, p4, d2, pd2)
    poly_mul1(u3, u4, p3, p4, d2, &pd2);
    //  call poly_sub1(d1, d2, pd1, pd2, u5, p5)
    poly_sub1(d1, d2, pd1, pd2, u5, p5);

    return;
    // end subroutine poly_mul_sub1
}
//!----------------------------------------------------------------------------
// subroutine poly_mul1(u1, u2, p1, p2, u3, p3)
void poly_mul1(double u1[17], double u2[17], int p1, int p2, double u3[17], int *p3) {
    //  implicit none
    //  real(dp), dimension(0:16), intent(in) :: u1, u2
    //  integer, intent(in) :: p1, p2
    //  real(dp), dimension(0:16), intent(out) :: u3
    //  integer, intent(out) :: p3
    //  integer :: i1, i2, i3
    int i, i1, i2, i3;
    //  real(dp) :: u1i
    double u1i;

    //  p3 = p1 + p2
    *p3 = p1 + p2;
    //  u3(:) = 0.0d0
    for (i = 0; i < 17; i++)
        u3[i] = 0.;

    //  do i1 = 0, p1
    for (i1 = 0; i1 <= p1; i1++) {
        //     u1i = u1(i1)
        u1i = u1[i1];
        //     do i2 = 0, p2
        for (i2 = 0; i2 <= p2; i2++) {
            //        i3 = i1 + i2
            i3 = i1 + i2;
            //        u3(i3) = u3(i3) + u1i*u2(i2)
            u3[i3] = u3[i3] + u1i * u2[i2];
            //     end do
        }
        //  end do
    }

    return;
    // end subroutine poly_mul1
}
//!----------------------------------------------------------------------------
// subroutine poly_sub1(u1, u2, p1, p2, u3, p3)
void poly_sub1(double u1[17], double u2[17], int p1, int p2, double u3[17], int *p3) {
    //  implicit none
    //  real(dp), dimension(0:16), intent(in) :: u1, u2
    //  integer, intent(in) :: p1, p2
    //  real(dp), dimension(0:16), intent(out) :: u3
    //  integer, intent(out) :: p3
    //  integer :: i
    int i;

    //  p3 = max(p1, p2)
    *p3 = max(p1, p2);

    //  do i = 0, p3
    for (i = 0; i <= *p3; i++) {
        //     if (i > p2) then
        //        u3(i) = u1(i)
        if (i > p2)
            u3[i] = u1[i];
        //     else if (i > p1) then
        //        u3(i) = -u2(i)
        else if (i > p1)
            u3[i] = -u2[i];
        //     else
        //        u3(i) = u1(i) - u2(i)
        else
            u3[i] = u1[i] - u2[i];
        //     end if
        //  end do
    }

    return;
    // end subroutine poly_sub1
}
//!----------------------------------------------------------------------------
// subroutine coord_from_poly_roots(n_soln, roots, r_n1, r_a1, r_a3, r_c3, r_soln_n, r_soln_a,
// r_soln_c)
void coord_from_poly_roots(int *n_soln, double roots[max_soln], const Vec3& r_n1, const Vec3& r_a1,
                           const Vec3& r_a3, const Vec3& r_c3,
                           std::array<Eigen::Matrix3d, max_soln>& r_soln_n,
                           std::array<Eigen::Matrix3d, max_soln> &r_soln_a,
                           std::array<Eigen::Matrix3d, max_soln>& r_soln_c) {
    Vec3 ex, ey, ez, b_a1a2, b_a3a2, r_tmp;
    std::array<Vec3, 3> p_s, s1, s2, p_t, t1, t2;
    std::array<Vec3, 3> p_s_c, s1_s, s2_s, p_t_c, t1_s, t2_s;
    double angle, sig1_init, half_tan[3];
    double cos_tau[4], sin_tau[4], cos_sig[3], sin_sig[3], ht, tmp, sig1;
    Vec3 r_s, r_t, r0;
    std::array<Vec3, 3> r_n, r_a, r_c;
    Eigen::Quaterniond p;
    Eigen::Matrix3d Us;
    int i_soln, i, j;
    double a1c1, c1n2, n2a2, a2c2, c2n3, n3a3, a1a2, a2a3;
    Vec3 rr_a1c1, rr_c1n2, rr_n2a2, rr_a2c2, rr_c2n3, rr_n3a3, rr_a1a2, rr_a2a3;
    double a3a1a2, a2a3a1, n1a1c1, n2a2c2, n3a3c3, a1c1n2a2, a2c2n3a3;
    double tmp_value;
    Vec3 ex_tmp;
    Vec3 tmp_array, tmp_array1, tmp_array2, tmp_array3;
    Vec3 mat1, mat2, mat3, mat4, mat5;
    Vec3 mat11, mat22, mat33, mat44, mat55;

    if (*n_soln == 0)
        return;
    for (i = 0; i < 3; i++)
        ex[i] = b_a1a3[i];
    cross(r_a1n1, ex, ez);
    tmp_value = sqrt(dot_product(ez, ez));
    for (i = 0; i < 3; i++)
        ez[i] = ez[i] / tmp_value;
    cross(ez, ex, ey);
    //  ! vertual bond vectors in the reference plane
    for (i = 0; i < 3; i++) {
        b_a1a2[i] = -cos_alpha[0] * ex[i] + sin_alpha[0] * ey[i];
        b_a3a2[i] = cos_alpha[2] * ex[i] + sin_alpha[2] * ey[i];
    }
    //  !! Define cone coordinates for each angle joint.
    //  ! (p_s,s1,s2) and (p_t,t1,t2):  Right Orthonormal systems
    p_s[0] -= ex;
    s1[0] = ez;
    s2[0] = ey;
    p_t[0] = b_a1a2;
    t1[0]  = ez;
    t2[0]  = sin_alpha[0] * ex + cos_alpha[0] * ey;

    //  ! residue 2
    p_s[1] = -b_a1a2;
    s1[1]  = -ez;
    s2[1]  = t2[0];
    p_t[1] = -b_a3a2;
    t1[1]  = -ez;
    t2[1]  = sin_alpha[2] * ex - cos_alpha[2] * ey;

    //  ! residue 3
    p_s[2] = b_a3a2;
    s1[2]  = ez;
    s2[2]  = t2[1];
    p_t[2] = ex;
    t1[2]  = ez;
    t2[2]  = -ey;

    //  ! scale vectors
    for (i = 0; i < 3; i++) {
        p_s_c[i] = p_s[i] * cos_xi[i];
        s1_s[i]  = s1[i] * sin_xi[i];
        s2_s[i]  = s2[i] * sin_xi[i];
        p_t_c[i] = p_t[i] * cos_eta[i];
        t1_s[i]  = t1[i] * sin_eta[i];
        t2_s[i]  = t2[i] * sin_eta[i];
    }

    //  ! initial sig(1)
    r_tmp = (r_a1n1 / len_na[0] - p_s_c[0]) / sin_xi[0];
    calc_bnd_ang(s1[0], r_tmp, &angle);
    sig1_init = sign(angle, dot_product(r_tmp, s2[0]));

    //  ! CA
    r_a[0] = r_a1;
    r_a[1] = r_a1 + len_aa[1] * b_a1a2;
    r_a[2] = r_a3;
    r0 = r_a1;

    for (i_soln = 0; i_soln < *n_soln; i_soln++) {
        half_tan[2] = roots[i_soln];
        half_tan[1] = calc_t2(half_tan[2]);
        half_tan[0] = calc_t1(half_tan[2], half_tan[1]);
        for (i = 1; i <= 3; i++) {
            ht         = half_tan[i - 1];
            tmp        = 1.0e0 + ht * ht;
            cos_tau[i] = (1.0e0 - ht * ht) / tmp;
            sin_tau[i] = 2.0e0 * ht / tmp;
        }
        cos_tau[0] = cos_tau[3];
        sin_tau[0] = sin_tau[3];
        for (i = 0; i < 3; i++) {
            cos_sig[i] = cos_delta[i] * cos_tau[i] + sin_delta[i] * sin_tau[i];
            sin_sig[i] = sin_delta[i] * cos_tau[i] - cos_delta[i] * sin_tau[i];
        }
        for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++) {
                r_s[j]    = p_s_c[i][j] + cos_sig[i] * s1_s[i][j] + sin_sig[i] * s2_s[i][j];
                r_t[j]    = p_t_c[i][j] + cos_tau[i + 1] * t1_s[i][j] + sin_tau[i + 1] * t2_s[i][j];
                r_n[i][j] = r_s[j] * len_na[i] + r_a[i][j];
                r_c[i][j] = r_t[j] * len_ac[i] + r_a[i][j];
            }

        //     ! rotate back atoms by -(sig(1) - sig1_init) around -ex
        sig1 = atan2(sin_sig[0], cos_sig[0]);
        ex_tmp -= ex;
        tmp_value = -(sig1 - sig1_init) * 0.25;
        quaternion(ex_tmp, tmp_value, p);
        rotation_matrix(p, Us);
        mat11 = r_c[0] - r0;
        mat22 = r_n[1] - r0;
        mat33 = r_a[1] - r0;
        mat44 = r_c[1] - r0;
        mat55 = r_n[2] - r0;
        matmul(Us, mat11, mat1);
        matmul(Us, mat22, mat2);
        matmul(Us, mat33, mat3);
        matmul(Us, mat44, mat4);
        matmul(Us, mat55, mat5);

        r_soln_n[i_soln].col(0) = r_n1;
        r_soln_a[i_soln].col(0) = r_a1;
        r_soln_c[i_soln].col(0) = mat1 + r0;
        r_soln_n[i_soln].col(1) = mat2 + r0;
        r_soln_a[i_soln].col(1) = mat3 + r0;
        r_soln_c[i_soln].col(1) = mat4 + r0;
        r_soln_n[i_soln].col(2) = mat5 + r0;
        r_soln_a[i_soln].col(2) = r_a3;
        r_soln_c[i_soln].col(2) = r_c3;
        //  end do
    }

    return;
    // end subroutine coord_from_poly_roots
}
//!-----------------------------------------------------------------------
// function calc_t2(t0)
double calc_t2(double t0) {
    //  implicit none
    //  real(dp), intent(in) :: t0
    //  real(dp) :: calc_t2
    double tmp_value;
    //  real(dp) :: B0, B1, B2, A0, A1, A2, A3, A4, B2_2, B2_3
    double B0, B1, B2, A0, A1, A2, A3, A4, B2_2, B2_3;
    //  real(dp) :: K0, K1, K2, K3, t0_2, t0_3, t0_4
    double K0, K1, K2, K3, t0_2, t0_3, t0_4;

    //  t0_2 = t0*t0
    t0_2 = t0 * t0;
    //  t0_3 = t0_2*t0
    t0_3 = t0_2 * t0;
    //  t0_4 = t0_3*t0
    t0_4 = t0_3 * t0;

    //  A0 = Q(0,0) + Q(1,0)*t0 + Q(2,0)*t0_2 + Q(3,0)*t0_3 + Q(4,0)*t0_4
    A0 = Q[0][0] + Q[0][1] * t0 + Q[0][2] * t0_2 + Q[0][3] * t0_3 + Q[0][4] * t0_4;
    //  A1 = Q(0,1) + Q(1,1)*t0 + Q(2,1)*t0_2 + Q(3,1)*t0_3 + Q(4,1)*t0_4
    A1 = Q[1][0] + Q[1][1] * t0 + Q[1][2] * t0_2 + Q[1][3] * t0_3 + Q[1][4] * t0_4;
    //  A2 = Q(0,2) + Q(1,2)*t0 + Q(2,2)*t0_2 + Q(3,2)*t0_3 + Q(4,2)*t0_4
    A2 = Q[2][0] + Q[2][1] * t0 + Q[2][2] * t0_2 + Q[2][3] * t0_3 + Q[2][4] * t0_4;
    //  A3 = Q(0,3) + Q(1,3)*t0 + Q(2,3)*t0_2 + Q(3,3)*t0_3 + Q(4,3)*t0_4
    A3 = Q[3][0] + Q[3][1] * t0 + Q[3][2] * t0_2 + Q[3][3] * t0_3 + Q[3][4] * t0_4;
    //  A4 = Q(0,4) + Q(1,4)*t0 + Q(2,4)*t0_2 + Q(3,4)*t0_3 + Q(4,4)*t0_4
    A4 = Q[4][0] + Q[4][1] * t0 + Q[4][2] * t0_2 + Q[4][3] * t0_3 + Q[4][4] * t0_4;

    //  B0 = R(0,0) + R(1,0)*t0 + R(2,0)*t0_2
    B0 = R[0][0] + R[0][1] * t0 + R[0][2] * t0_2;
    //  B1 = R(0,1) + R(1,1)*t0 + R(2,1)*t0_2
    B1 = R[1][0] + R[1][1] * t0 + R[1][2] * t0_2;
    //  B2 = R(0,2) + R(1,2)*t0 + R(2,2)*t0_2
    B2 = R[2][0] + R[2][1] * t0 + R[2][2] * t0_2;

    //  B2_2 = B2*B2
    B2_2 = B2 * B2;
    //  B2_3 = B2_2*B2
    B2_3 = B2_2 * B2;

    //  K0 = A2*B2 - A4*B0
    K0 = A2 * B2 - A4 * B0;
    //  K1 = A3*B2 - A4*B1
    K1 = A3 * B2 - A4 * B1;
    //  K2 = A1*B2_2 - K1*B0
    K2 = A1 * B2_2 - K1 * B0;
    //  K3 = K0*B2 - K1*B1
    K3 = K0 * B2 - K1 * B1;

    //  calc_t2 = (K3*B0 - A0*B2_3)/(K2*B2 - K3*B1)
    tmp_value = (K3 * B0 - A0 * B2_3) / (K2 * B2 - K3 * B1);

    return tmp_value;
    // end function calc_t2
}
//!-----------------------------------------------------------------------
// function calc_t1(t0, t2)
double calc_t1(double t0, double t2) {
    //  implicit none
    //  real(dp), intent(in) :: t0, t2
    //  real(dp) :: calc_t1
    double tmp_value;
    //  real(dp) :: U11, U12, U13, U31, U32, U33
    double U11, U12, U13, U31, U32, U33;
    //  real(dp) :: t0_2, t2_2
    double t0_2, t2_2;

    //  t0_2 = t0*t0
    t0_2 = t0 * t0;
    //  t2_2 = t2*t2
    t2_2 = t2 * t2;

    //  U11 = C0(0,1) + C0(1,1)*t0 + C0(2,1)*t0_2
    U11 = C0[0][0] + C0[0][1] * t0 + C0[0][2] * t0_2;
    //  U12 = C1(0,1) + C1(1,1)*t0 + C1(2,1)*t0_2
    U12 = C1[0][0] + C1[0][1] * t0 + C1[0][2] * t0_2;
    //  U13 = C2(0,1) + C2(1,1)*t0 + C2(2,1)*t0_2
    U13 = C2[0][0] + C2[0][1] * t0 + C2[0][2] * t0_2;
    //  U31 = C0(0,2) + C0(1,2)*t2 + C0(2,2)*t2_2
    U31 = C0[1][0] + C0[1][1] * t2 + C0[1][2] * t2_2;
    //  U32 = C1(0,2) + C1(1,2)*t2 + C1(2,2)*t2_2
    U32 = C1[1][0] + C1[1][1] * t2 + C1[1][2] * t2_2;
    //  U33 = C2(0,2) + C2(1,2)*t2 + C2(2,2)*t2_2
    U33 = C2[1][0] + C2[1][1] * t2 + C2[1][2] * t2_2;

    //  calc_t1 = (U31*U13-U11*U33)/(U12*U33-U13*U32)
    tmp_value = (U31 * U13 - U11 * U33) / (U12 * U33 - U13 * U32);

    return tmp_value;
    // end function calc_t1
}

void calc_dih_ang(const Vec3& r1, const Vec3& r2, const Vec3& r3, double *angle) {
    // Compute plane normals
    Vec3 p, q;
    CrossProduct(r1, r2, &p); // Normal to plane abc
    CrossProduct(r2, r3, &q); // Normal to plane bcd

    // Calculate the cosine of the angle (dot product of normals)
    double cos_phi = Dot(p, q) / sqrt(Dot(p, p) * Dot(q, q));
    
    // To get the sign and use atan2
    Vec3 s;
    CrossProduct(p, r2, &s);
    double sin_phi = Dot(s, q) / (sqrt(Dot(s, s) * Dot(q, q)));

    // atan2 is stable and handles the quadrants automatically
    *angle = atan2(sin_phi, cos_phi);
}

double struct_calc_dih_ang(const Vec3& r1, const Vec3& r2, const Vec3& r3) {
    // NOTE: This is duplicate of calc_dih_ang but returns a float for use in the struct version of the code.
    // Compute plane normals
    Vec3 p, q;
    CrossProduct(r1, r2, &p); // Normal to plane abc
    CrossProduct(r2, r3, &q); // Normal to plane bcd

    // Calculate the cosine of the angle (dot product of normals)
    double cos_phi = Dot(p, q) / sqrt(Dot(p, p) * Dot(q, q));
    
    // To get the sign and use atan2
    Vec3 s;
    CrossProduct(p, r2, &s);
    double sin_phi = Dot(s, q) / (sqrt(Dot(s, s) * Dot(q, q)));

    return atan2(sin_phi, cos_phi);
}
//!-----------------------------------------------------------------------
// subroutine calc_bnd_ang(r1, r2, angle)
void calc_bnd_ang(const Vec3& r1, const Vec3& r2, double *angle) {
    double arg = r1.dot(r2);
    arg = std::clamp(arg, -1.0, 1.0); // for numerical safety
    *angle = std::acos(arg);
    return;
}

void c_bnd_len(const Vec3& r1, const Vec3& r2, double *length) {
    *length = (r1 - r2).norm();
    return;
}

void c_bnd_ang(const Vec3& r1, const Vec3& r2, const Vec3& r3, double *angle) {
    Vec3 r21 = (r1 - r2).normalized();
    Vec3 r23 = (r3 - r2).normalized();
    calc_bnd_ang(r21, r23, angle);
}

void c_dih_ang(const Vec3& r1, const Vec3& r2, const Vec3& r3, const Vec3& r4, double *angle) {
    Vec3 r12, r23, r34;
    calc_dih_ang(r12, r23, r34, angle);
    return;
}
//!-----------------------------------------------------------------------
// subroutine cross(p, q, s)
void cross(const Vec3& p, const Vec3& q, Vec3& s) {
    CrossProduct(p, q, &s);
    return;
}
//!-----------------------------------------------------------------------
// subroutine quaternion(axis, quarter_ang, p)
void quaternion(const Vec3& axis, double quarter_ang, Eigen::Quaterniond& p) {
    double rotation_angle = 4.0 * quarter_ang;
    Eigen::AngleAxisd aa(rotation_angle, axis.normalized());
    p = aa;
}

//!-----------------------------------------------------------------------
// subroutine rotation_matrix(q, U)
void rotation_matrix(const Eigen::Quaterniond& q, Eigen::Matrix3d& U) {
    U = q.toRotationMatrix();
    return;
}
//!----------------------------------------------------------------------------
// END MODULE tripep_closure
//!----------------------------------------------------------------------------
void apply_rotation_eigen(const Eigen::Ref<Vec3> obj_i, 
                          Eigen::Ref<Vec3> obj_f, 
                          const Eigen::Ref<Vec3> tor2, 
                          const Eigen::Ref<Vec3> tor3, 
                          double total_angle) {
    // Normalized rotation axis (tor3 - tor2)
    Vec3 axis = (tor3 - tor2).normalized();
    // Create the rotation object (AngleAxis)
    Eigen::AngleAxisd rotation(total_angle, axis);
    // Apply rotation: Translate to origin (relative to tor3), rotate, translate back
    obj_f = rotation * (obj_i - tor3) + tor3;
}


void get_rot(const Eigen::Ref<Vec3> obj_i, Eigen::Ref<Vec3> obj_f, 
             const Eigen::Ref<Vec3> tor1, const Eigen::Ref<Vec3> tor2, const Eigen::Ref<Vec3> tor3,
             const Eigen::Ref<Vec3> tor4_i, const Eigen::Ref<Vec3> tor4_f) {
    double ang_i, ang_f;
    c_dih_ang(tor1, tor2, tor3, tor4_i, &ang_i);
    c_dih_ang(tor1, tor2, tor3, tor4_f, &ang_f);
    double total_angle = ang_f - ang_i;
    apply_rotation_eigen(obj_i, obj_f, tor2, tor3, total_angle);
}
void driver_rot(const Eigen::Ref<Vec3> obj_i, Eigen::Ref<Vec3> obj_f, 
                const Eigen::Ref<Vec3> tor2, const Eigen::Ref<Vec3> tor3,
                double dih_change) {
    apply_rotation_eigen(obj_i, obj_f, tor2, tor3, dih_change);
}