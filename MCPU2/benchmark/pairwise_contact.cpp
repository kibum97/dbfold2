#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <chrono>

struct Atom {
    Eigen::Map<Eigen::Vector3d> position;
    Atom(double* raw_ptr) : position(raw_ptr) {}
};

double compute_distance_sum(const std::vector<Atom>& atoms) {
    double sum = 0.0;
    for (size_t i = 0; i < atoms.size(); ++i) {
        for (size_t j = i + 1; j < atoms.size(); ++j) {
            sum += (atoms[i].position - atoms[j].position).norm();
        }
    }
    return sum;
}

double compute_distance_sum_direct(const Eigen::MatrixXd& positions) {
    double sum = 0.0;
    int N = positions.cols();
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            sum += (positions.col(i) - positions.col(j)).norm();
        }
    }
    return sum;
}

int main() {
    constexpr int N = 1000;  // number of atoms

    Eigen::MatrixXd positions = Eigen::MatrixXd::Random(3, N);

    // Create atoms with mapped columns
    std::vector<Atom> atoms;
    atoms.reserve(N);
    for (int i = 0; i < N; ++i) {
        atoms.emplace_back(positions.col(i).data());
    }

    // Benchmark 1: Atom-based access
    auto start_atom = std::chrono::high_resolution_clock::now();
    double sum1 = compute_distance_sum(atoms);
    auto end_atom = std::chrono::high_resolution_clock::now();

    // Benchmark 2: Direct access
    auto start_direct = std::chrono::high_resolution_clock::now();
    double sum2 = compute_distance_sum_direct(positions);
    auto end_direct = std::chrono::high_resolution_clock::now();

    // Results
    auto duration_atom = std::chrono::duration_cast<std::chrono::milliseconds>(end_atom - start_atom).count();
    auto duration_direct = std::chrono::duration_cast<std::chrono::milliseconds>(end_direct - start_direct).count();

    std::cout << "Atom struct Map access time: " << duration_atom << " ms, Sum: " << sum1 << std::endl;
    std::cout << "Direct matrix col access time: " << duration_direct << " ms, Sum: " << sum2 << std::endl;

    return 0;
}
