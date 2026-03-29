#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <chrono>

int main() {
    constexpr int N = 3000; // number of BB atoms
    constexpr int M = 5000; // number of SC atoms

    // Random seed for reproducibility
    srand(42);

    // Initialize matrices
    Eigen::Matrix3Xf BB = Eigen::Matrix3Xf::Random(3, N);
    Eigen::Matrix3Xf SC = Eigen::Matrix3Xf::Random(3, M);

    // Concatenate BB and SC into one matrix: all = [BB, SC]
    Eigen::Matrix3Xf all(3, N + M);
    all << BB, SC;

    // ===========================
    // Vectorized version
    // ===========================
    auto t1 = std::chrono::high_resolution_clock::now();

    Eigen::RowVectorXf all_sq = all.colwise().squaredNorm();
    Eigen::RowVectorXf sc_sq = SC.colwise().squaredNorm();

    Eigen::MatrixXf dist2 = all_sq.transpose().replicate(1, M)
                          + sc_sq.replicate(N + M, 1)
                          - 2 * (all.transpose() * SC);

    Eigen::MatrixXf dist_matrix = dist2.cwiseMax(0).cwiseSqrt();

    auto t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration_vectorized = t2 - t1;
    std::cout << "Vectorized time: " << duration_vectorized.count() << " s\n";


    // ===========================
    // Loop version
    // ===========================
    auto t3 = std::chrono::high_resolution_clock::now();

    Eigen::MatrixXf loop_matrix(N + M, M);

    for (int i = 0; i < N + M; ++i) {
        for (int j = 0; j < M; ++j) {
            float dist = (all.col(i) - SC.col(j)).norm();
            loop_matrix(i, j) = dist;
        }
    }

    auto t4 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration_loop = t4 - t3;
    std::cout << "Looped time: " << duration_loop.count() << " s\n";


    // Optional: check that the matrices are numerically similar
    float max_diff = (dist_matrix - loop_matrix).cwiseAbs().maxCoeff();
    std::cout << "Max difference: " << max_diff << std::endl;

    return 0;
}