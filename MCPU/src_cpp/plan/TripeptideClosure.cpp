#include "TripeptideClosure.hpp"
#include <iostream>

namespace TLC {

// High-performance replacement for calc_dih_ang
double TripeptideSolver::calc_dih_ang(const Eigen::Vector3d& r1, const Eigen::Vector3d& r2, const Eigen::Vector3d& r3) {
    Eigen::Vector3d p = r1.cross(r2);
    Eigen::Vector3d q = r2.cross(r3);
    Eigen::Vector3d s = r3.cross(r1);
    
    double arg = p.dot(q) / std::sqrt(p.squaredNorm() * q.squaredNorm());
    arg = std::clamp(arg, -1.0, 1.0); // Safe bounding
    
    double angle = std::acos(arg);
    return (s.dot(r2) >= 0.0) ? angle : -angle;
}

// High-performance replacement for calc_bnd_ang
double TripeptideSolver::calc_bnd_ang(const Eigen::Vector3d& r1, const Eigen::Vector3d& r2) {
    double arg = r1.dot(r2); // Assumes r1, r2 are normalized
    arg = std::clamp(arg, -1.0, 1.0);
    return std::acos(arg);
}

void TripeptideSolver::initialize(const std::array<double, 6>& b_len, 
                                  const std::array<double, 7>& b_ang, 
                                  const std::array<double, 2>& t_ang) {
    len0 = b_len;
    b_ang0 = b_ang;
    t_ang0 = t_ang;

    Eigen::Vector3d axis(1.0, 0.0, 0.0);
    Eigen::Vector3d rr_c1 = Eigen::Vector3d::Zero();

    for (int i = 0; i < 2; i++) {
        Eigen::Vector3d rr_a1(std::cos(b_ang0[3*i+1]) * len0[3*i], 
                              std::sin(b_ang0[3*i+1]) * len0[3*i], 0.0);
                              
        Eigen::Vector3d rr_n2(len0[3*i+1], 0.0, 0.0);
        Eigen::Vector3d rr_c1a1 = rr_a1 - rr_c1;
        
        Eigen::Vector3d rr_n2a2_ref(-std::cos(b_ang0[3*i+2]) * len0[3*i+2], 
                                     std::sin(b_ang0[3*i+2]) * len0[3*i+2], 0.0);

        // Replace quaternion & rotation_matrix with Eigen::AngleAxisd
        // Old code used angle * 0.25, which corresponds to rotation of angle/2.0
        Eigen::Matrix3d Us = Eigen::AngleAxisd(t_ang0[i] / 2.0, axis).toRotationMatrix();

        Eigen::Vector3d rr_a2 = (Us * rr_n2a2_ref) + rr_n2;
        Eigen::Vector3d rr_a1a2 = rr_a2 - rr_a1;
        
        double len1 = rr_a1a2.norm();
        len_aa[i+1] = len1;

        Eigen::Vector3d bb_c1a1 = rr_c1a1.normalized();
        Eigen::Vector3d bb_a1a2 = rr_a1a2.normalized();
        Eigen::Vector3d bb_a2n2 = (rr_n2 - rr_a2).normalized();

        xi[i+1] = calc_bnd_ang(-bb_a1a2, bb_a2n2);
        eta[i] = calc_bnd_ang(bb_a1a2, -bb_c1a1);
        
        delta[i+1] = PI - calc_dih_ang(bb_c1a1, bb_a1a2, bb_a2n2);
    }

    double a_min = b_ang[3] - (xi[1] + eta[1]);
    double a_max = std::min(b_ang[3] + (xi[1] + eta[1]), PI);

    aa13_min_sqr = std::pow(len_aa[1], 2) + std::pow(len_aa[2], 2) - 2.0 * len_aa[1] * len_aa[2] * std::cos(a_min);
    aa13_max_sqr = std::pow(len_aa[1], 2) + std::pow(len_aa[2], 2) - 2.0 * len_aa[1] * len_aa[2] * std::cos(a_max);
}

std::vector<Solution> TripeptideSolver::solve(const Eigen::Vector3d& r_n1, 
                                              const Eigen::Vector3d& r_a1, 
                                              const Eigen::Vector3d& r_a3, 
                                              const Eigen::Vector3d& r_c3) {
    std::vector<Solution> valid_solutions;
    
    Eigen::Vector3d b_a1a3, b_a1n1, b_a3c3;
    if (!get_input_angles(r_n1, r_a1, r_a3, r_c3, b_a1a3, b_a1n1, b_a3c3)) {
        return valid_solutions; // No geometric solution possible
    }

    // Stack allocated fixed size vector for speed
    Eigen::Matrix<double, 17, 1> poly_coeff;
    get_poly_coeff(poly_coeff);

    // ====================================================================
    // REPLACING STURM SOLVER WITH EIGEN POLYNOMIAL SOLVER
    // ====================================================================
    Eigen::PolynomialSolver<double, 16> solver;
    solver.compute(poly_coeff); // poly_coeff MUST be ordered from a0 to a16

    std::vector<double> real_roots;
    real_roots.reserve(16);
    
    for (int i = 0; i < solver.roots().size(); ++i) {
        // Filter out imaginary roots (tolerance 1e-9)
        if (std::abs(solver.roots()[i].imag()) < 1e-9) {
            real_roots.push_back(solver.roots()[i].real());
        }
    }

    if (real_roots.empty()) {
        return valid_solutions;
    }

    coord_from_poly_roots(real_roots, r_n1, r_a1, r_a3, r_c3, b_a1a3, valid_solutions);
    return valid_solutions;
}

bool TripeptideSolver::get_input_angles(const Eigen::Vector3d& r_n1, const Eigen::Vector3d& r_a1, 
                                        const Eigen::Vector3d& r_a3, const Eigen::Vector3d& r_c3,
                                        Eigen::Vector3d& b_a1a3, Eigen::Vector3d& b_a1n1, Eigen::Vector3d& b_a3c3) {
    
    Eigen::Vector3d r_a1a3 = r_a3 - r_a1;
    double dr_sqr = r_a1a3.squaredNorm();
    len_aa[0] = std::sqrt(dr_sqr);

    // Boundary check: If the anchors are too far or too close, no loop can bridge them.
    if (dr_sqr < aa13_min_sqr || dr_sqr > aa13_max_sqr) {
        return false; 
    }

    // Bond lengths
    Eigen::Vector3d r_a1n1 = r_n1 - r_a1;
    len_na[0] = r_a1n1.norm();
    len_na[1] = len0[2];
    len_na[2] = len0[5];

    Eigen::Vector3d r_a3c3 = r_c3 - r_a3;
    len_ac[0] = len0[0];
    len_ac[1] = len0[3];
    len_ac[2] = r_a3c3.norm();

    // Unit vectors
    b_a1n1 = r_a1n1 / len_na[0];
    b_a3c3 = r_a3c3 / len_ac[2];
    b_a1a3 = r_a1a3 / len_aa[0];

    // Calculate internal angles using the robust helper
    delta[3] = calc_dih_ang(-b_a1n1, b_a1a3, b_a3c3);
    delta[0] = delta[3];

    xi[0] = calc_bnd_ang(-b_a1a3, b_a1n1);
    eta[2] = calc_bnd_ang(b_a1a3, b_a3c3);

    // Precompute sines and cosines to avoid re-evaluating in the tight loop
    for (int i = 0; i < 3; i++) {
        cos_delta[i+1] = std::cos(delta[i+1]);
        sin_delta[i+1] = std::sin(delta[i+1]);
        
        cos_xi[i] = std::cos(xi[i]);
        sin_xi[i] = std::sin(xi[i]);
        
        cos_eta[i] = std::cos(eta[i]);
        sin_eta[i] = std::sin(eta[i]);
    }
    
    cos_delta[0] = cos_delta[3];
    sin_delta[0] = sin_delta[3];

    theta[0] = b_ang0[0];
    theta[1] = b_ang0[3];
    theta[2] = b_ang0[6];

    for (int i = 0; i < 3; i++) {
        cos_theta[i] = std::cos(theta[i]);
    }

    // Alpha angles (Law of Cosines)
    cos_alpha[0] = -(std::pow(len_aa[0], 2) + std::pow(len_aa[1], 2) - std::pow(len_aa[2], 2)) / (2.0 * len_aa[0] * len_aa[1]);
    alpha[0] = std::acos(cos_alpha[0]);
    sin_alpha[0] = std::sin(alpha[0]);

    cos_alpha[1] = (std::pow(len_aa[1], 2) + std::pow(len_aa[2], 2) - std::pow(len_aa[0], 2)) / (2.0 * len_aa[1] * len_aa[2]);
    alpha[1] = std::acos(cos_alpha[1]);
    sin_alpha[1] = std::sin(alpha[1]);

    alpha[2] = PI - alpha[0] + alpha[1];
    cos_alpha[2] = std::cos(alpha[2]);
    sin_alpha[2] = std::sin(alpha[2]);

    return true; 
}

void TripeptideSolver::coord_from_poly_roots(const std::vector<double>& roots, 
                                             const Eigen::Vector3d& r_n1, const Eigen::Vector3d& r_a1, 
                                             const Eigen::Vector3d& r_a3, const Eigen::Vector3d& r_c3,
                                             const Eigen::Vector3d& b_a1a3,
                                             std::vector<Solution>& solutions) {
    // Define body frame (ex, ey, ez)
    Eigen::Vector3d ex = b_a1a3;
    Eigen::Vector3d r_a1n1 = r_n1 - r_a1;
    Eigen::Vector3d ez = r_a1n1.cross(ex).normalized();
    Eigen::Vector3d ey = ez.cross(ex);

    // Virtual bond vectors in reference plane
    Eigen::Vector3d b_a1a2 = -cos_alpha[0]*ex + sin_alpha[0]*ey;
    Eigen::Vector3d b_a3a2 = cos_alpha[2]*ex + sin_alpha[2]*ey;

    // Set up local cone coordinate frames (p_s, s1, s2) and (p_t, t1, t2)
    std::array<Eigen::Vector3d, 3> p_s, s1, s2, p_t, t1, t2;
    std::array<Eigen::Vector3d, 3> p_s_c, s1_s, s2_s, p_t_c, t1_s, t2_s;

    // Residue 1
    p_s[0] = -ex;       s1[0] = ez;   s2[0] = ey;
    p_t[0] = b_a1a2;    t1[0] = ez;   t2[0] = sin_alpha[0]*ex + cos_alpha[0]*ey;

    // Residue 2
    p_s[1] = -b_a1a2;   s1[1] = -ez;  s2[1] = t2[0];
    p_t[1] = -b_a3a2;   t1[1] = -ez;  t2[1] = sin_alpha[2]*ex - cos_alpha[2]*ey;

    // Residue 3
    p_s[2] = b_a3a2;    s2[2] = t2[1]; s1[2] = ez;
    p_t[2] = ex;        t1[2] = ez;    t2[2] = -ey;

    // Pre-scale vectors by xi and eta to save flops in the inner loop
    for(int i=0; i<3; ++i) {
        p_s_c[i] = p_s[i] * cos_xi[i];
        s1_s[i]  = s1[i]  * sin_xi[i];
        s2_s[i]  = s2[i]  * sin_xi[i];
        p_t_c[i] = p_t[i] * cos_eta[i];
        t1_s[i]  = t1[i]  * sin_eta[i];
        t2_s[i]  = t2[i]  * sin_eta[i];
    }

    // Initial sig(1)
    Eigen::Vector3d r_tmp = (r_a1n1 / len_na[0] - p_s_c[0]) / sin_xi[0];
    double angle = calc_bnd_ang(s1[0], r_tmp);
    double sig1_init = (r_tmp.dot(s2[0]) >= 0.0) ? std::abs(angle) : -std::abs(angle);

    std::array<Eigen::Vector3d, 3> r_a;
    r_a[0] = r_a1;
    r_a[1] = r_a1 + len_aa[1] * b_a1a2;
    r_a[2] = r_a3;
    Eigen::Vector3d r0 = r_a1;

    for (double root : roots) {
        Solution sol;
        std::array<double, 3> half_tan;
        half_tan[2] = root;
        half_tan[1] = calc_t2(half_tan[2]);
        half_tan[0] = calc_t1(half_tan[2], half_tan[1]);

        std::array<double, 4> cos_tau, sin_tau;
        std::array<double, 3> cos_sig, sin_sig;

        for (int i = 1; i <= 3; i++) {
            double ht = half_tan[i-1];
            double tmp = 1.0 + ht*ht;
            cos_tau[i] = (1.0 - ht*ht) / tmp;
            sin_tau[i] = 2.0 * ht / tmp;
        }
        cos_tau[0] = cos_tau[3];
        sin_tau[0] = sin_tau[3];

        for (int i = 0; i < 3; i++) {
            cos_sig[i] = cos_delta[i]*cos_tau[i] + sin_delta[i]*sin_tau[i];
            sin_sig[i] = sin_delta[i]*cos_tau[i] - cos_delta[i]*sin_tau[i];
        }

        std::array<Eigen::Vector3d, 3> r_n, r_c;
        for (int i = 0; i < 3; i++) {
            Eigen::Vector3d r_s = p_s_c[i] + cos_sig[i]*s1_s[i] + sin_sig[i]*s2_s[i];
            Eigen::Vector3d r_t = p_t_c[i] + cos_tau[i+1]*t1_s[i] + sin_tau[i+1]*t2_s[i];
            
            r_n[i] = r_s * len_na[i] + r_a[i];
            r_c[i] = r_t * len_ac[i] + r_a[i];
        }

        // Rotate back atoms around -ex
        double sig1 = std::atan2(sin_sig[0], cos_sig[0]);
        
        // Eigen::AngleAxisd expects (angle, axis). 
        // Original code used quaternion(axis, angle*0.25), which is a half-angle transformation.
        // The equivalent standard rotation is `-(sig1 - sig1_init) / 2.0`.
        Eigen::AngleAxisd rotation(-(sig1 - sig1_init) / 2.0, -ex);
        Eigen::Matrix3d Us = rotation.toRotationMatrix();

        sol.r_n[0] = r_n1;
        sol.r_a[0] = r_a1;
        sol.r_c[0] = (Us * (r_c[0] - r0)) + r0;
        
        sol.r_n[1] = (Us * (r_n[1] - r0)) + r0;
        sol.r_a[1] = (Us * (r_a[1] - r0)) + r0;
        sol.r_c[1] = (Us * (r_c[1] - r0)) + r0;
        
        sol.r_n[2] = (Us * (r_n[2] - r0)) + r0;
        sol.r_a[2] = r_a3;
        sol.r_c[2] = r_c3;

        solutions.push_back(sol);
    }
}

double TripeptideSolver::calc_t2(double t0) const {
    double t0_2 = t0 * t0;
    double t0_3 = t0_2 * t0;
    double t0_4 = t0_3 * t0;

    double A0 = Q(0,0) + Q(1,0)*t0 + Q(2,0)*t0_2 + Q(3,0)*t0_3 + Q(4,0)*t0_4;
    double A1 = Q(0,1) + Q(1,1)*t0 + Q(2,1)*t0_2 + Q(3,1)*t0_3 + Q(4,1)*t0_4;
    double A2 = Q(0,2) + Q(1,2)*t0 + Q(2,2)*t0_2 + Q(3,2)*t0_3 + Q(4,2)*t0_4;
    double A3 = Q(0,3) + Q(1,3)*t0 + Q(2,3)*t0_2 + Q(3,3)*t0_3 + Q(4,3)*t0_4;
    double A4 = Q(0,4) + Q(1,4)*t0 + Q(2,4)*t0_2 + Q(3,4)*t0_3 + Q(4,4)*t0_4;

    double B0 = R(0,0) + R(1,0)*t0 + R(2,0)*t0_2;
    double B1 = R(0,1) + R(1,1)*t0 + R(2,1)*t0_2;
    double B2 = R(0,2) + R(1,2)*t0 + R(2,2)*t0_2;

    double B2_2 = B2 * B2;
    double B2_3 = B2_2 * B2;

    double K0 = A2*B2 - A4*B0;
    double K1 = A3*B2 - A4*B1;
    double K2 = A1*B2_2 - K1*B0;
    double K3 = K0*B2 - K1*B1;

    return (K3*B0 - A0*B2_3) / (K2*B2 - K3*B1);
}

double TripeptideSolver::calc_t1(double t0, double t2) const {
    double t0_2 = t0 * t0;
    double t2_2 = t2 * t2;

    double U11 = C0(0,1) + C1(0,1)*t0 + C2(0,1)*t0_2; // Adjusted index matching C0 structure
    double U12 = C0(1,1) + C1(1,1)*t0 + C2(1,1)*t0_2;
    double U13 = C0(2,1) + C1(2,1)*t0 + C2(2,1)*t0_2;
    
    double U31 = C0(0,2) + C1(0,2)*t2 + C2(0,2)*t2_2;
    double U32 = C0(1,2) + C1(1,2)*t2 + C2(1,2)*t2_2;
    double U33 = C0(2,2) + C1(2,2)*t2 + C2(2,2)*t2_2;

    return (U31*U13 - U11*U33) / (U12*U33 - U13*U32);
}

// Polynomial Multiplication: u3 = u1 * u2
void TripeptideSolver::poly_mul1(const Eigen::VectorXd& u1, const Eigen::VectorXd& u2, Eigen::VectorXd& u3) {
    int deg1 = u1.size() - 1;
    int deg2 = u2.size() - 1;
    
    u3 = Eigen::VectorXd::Zero(deg1 + deg2 + 1);
    
    for (int i = 0; i <= deg1; ++i) {
        for (int j = 0; j <= deg2; ++j) {
            u3[i + j] += u1[i] * u2[j];
        }
    }
}

// Polynomial Subtraction: u3 = u1 - u2
void TripeptideSolver::poly_sub1(const Eigen::VectorXd& u1, const Eigen::VectorXd& u2, Eigen::VectorXd& u3) {
    int deg1 = u1.size() - 1;
    int deg2 = u2.size() - 1;
    int max_deg = std::max(deg1, deg2);
    
    u3 = Eigen::VectorXd::Zero(max_deg + 1);
    
    for (int i = 0; i <= max_deg; ++i) {
        double val1 = (i <= deg1) ? u1[i] : 0.0;
        double val2 = (i <= deg2) ? u2[i] : 0.0;
        u3[i] = val1 - val2;
    }
}

// Combined Multiply and Subtract: u5 = (u1 * u2) - (u3 * u4)
// This is heavily used in the determinant expansion in `get_poly_coeff`
void TripeptideSolver::poly_mul_sub1(const Eigen::VectorXd& u1, const Eigen::VectorXd& u2, 
                                     const Eigen::VectorXd& u3, const Eigen::VectorXd& u4, 
                                     Eigen::VectorXd& u5) {
    Eigen::VectorXd d1, d2;
    poly_mul1(u1, u2, d1);
    poly_mul1(u3, u4, d2);
    poly_sub1(d1, d2, u5);
}

}