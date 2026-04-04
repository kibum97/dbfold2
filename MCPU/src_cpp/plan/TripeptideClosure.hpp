#pragma once

#include <Eigen/Dense>
#include <unsupported/Eigen/Polynomials>
#include <array>
#include <vector>
#include <cmath>
#include <algorithm>

namespace TLC {

constexpr int MAX_SOLN = 16;
constexpr int DEG_POL = 16;
constexpr double PI = 3.14159265358979323846;

// Structure to hold a single loop closure solution
struct Solution {
    std::array<Eigen::Vector3d, 3> r_n; // N atoms for residue 1, 2, 3
    std::array<Eigen::Vector3d, 3> r_a; // CA atoms for residue 1, 2, 3
    std::array<Eigen::Vector3d, 3> r_c; // C atoms for residue 1, 2, 3
};

class TripeptideSolver {
public:
    TripeptideSolver() = default;

    // Equivalent to initialize_loop_closure
    void initialize(const std::array<double, 6>& b_len, 
                    const std::array<double, 7>& b_ang, 
                    const std::array<double, 2>& t_ang);

    // Equivalent to solve_3pep_poly
    // Returns a vector of valid 3D conformations
    std::vector<Solution> solve(const Eigen::Vector3d& r_n1, 
                                const Eigen::Vector3d& r_a1, 
                                const Eigen::Vector3d& r_a3, 
                                const Eigen::Vector3d& r_c3);

private:
    // Internal State (Previously Globals)
    std::array<double, 6> len0;
    std::array<double, 7> b_ang0;
    std::array<double, 2> t_ang0;

    double aa13_min_sqr, aa13_max_sqr;
    std::array<double, 4> delta;
    std::array<double, 3> xi, eta, alpha, theta;
    std::array<double, 3> cos_alpha, sin_alpha, cos_theta, sin_theta;
    std::array<double, 4> cos_delta, sin_delta;
    std::array<double, 3> cos_xi, cos_eta, sin_xi, sin_eta;
    
    std::array<double, 3> len_na, len_ac, len_aa;

    // Polynomial coefficients matrix (fixed size, stack allocated)
    Eigen::Matrix<double, 3, 3> C0, C1, C2;
    Eigen::Matrix<double, 5, 17> Q;
    Eigen::Matrix<double, 3, 17> R;

    // Core subroutines translated
    bool get_input_angles(const Eigen::Vector3d& r_n1, const Eigen::Vector3d& r_a1, 
                          const Eigen::Vector3d& r_a3, const Eigen::Vector3d& r_c3,
                          Eigen::Vector3d& b_a1a3, Eigen::Vector3d& b_a1n1, Eigen::Vector3d& b_a3c3);
                          
    void get_poly_coeff(Eigen::Matrix<double, 17, 1>& poly_coeff);
    
    void coord_from_poly_roots(const std::vector<double>& roots, 
                               const Eigen::Vector3d& r_n1, const Eigen::Vector3d& r_a1, 
                               const Eigen::Vector3d& r_a3, const Eigen::Vector3d& r_c3,
                               const Eigen::Vector3d& b_a1a3,
                               std::vector<Solution>& solutions);

    // Helper Math Functions utilizing Eigen
    double calc_t2(double t0) const;
    double calc_t1(double t0, double t2) const;
    
    static double calc_dih_ang(const Eigen::Vector3d& r1, const Eigen::Vector3d& r2, const Eigen::Vector3d& r3);
    static double calc_bnd_ang(const Eigen::Vector3d& r1, const Eigen::Vector3d& r2);

    // Polynomial math helpers (using Eigen vectors for fixed size ops)
    void poly_mul1(const Eigen::VectorXd& u1, const Eigen::VectorXd& u2, Eigen::VectorXd& u3);
    void poly_sub1(const Eigen::VectorXd& u1, const Eigen::VectorXd& u2, Eigen::VectorXd& u3);
    void poly_mul_sub1(const Eigen::VectorXd& u1, const Eigen::VectorXd& u2, 
                       const Eigen::VectorXd& u3, const Eigen::VectorXd& u4, 
                       Eigen::VectorXd& u5);
};

} // namespace TLC