#include <iostream>
#include <vector>
#include <chrono>
#include <Eigen/Dense>

constexpr int ROWS = 1000;
constexpr int COLS = 1000;
constexpr int NUM_REPEATS = 10;

float benchmark_eigen() {
    Eigen::MatrixXf mat = Eigen::MatrixXf::Random(ROWS, COLS);
    float sum = 0.0f;
    for (int r = 0; r < NUM_REPEATS; ++r) {
        sum = 0.0f;
        for (int j = 0; j < COLS; ++j)
            for (int i = 0; i < ROWS; ++i)
                sum += mat(i, j);  // Column-major access
    }
    return sum;
}

float benchmark_vector_of_vector() {
    std::vector<std::vector<float>> mat(ROWS, std::vector<float>(COLS, 1.0f));
    float sum = 0.0f;
    for (int r = 0; r < NUM_REPEATS; ++r) {
        sum = 0.0f;
        for (int j = 0; j < COLS; ++j)
            for (int i = 0; i < ROWS; ++i)
                sum += mat[i][j];
    }
    return sum;
}

float benchmark_flattened_vector() {
    std::vector<float> mat(ROWS * COLS, 1.0f);
    float sum = 0.0f;
    for (int r = 0; r < NUM_REPEATS; ++r) {
        sum = 0.0f;
        for (int j = 0; j < COLS; ++j)
            for (int i = 0; i < ROWS; ++i)
                sum += mat[i * COLS + j];  // Row-major order
    }
    return sum;
}

int main() {
    using namespace std::chrono;

    auto t1 = high_resolution_clock::now();
    float s1 = benchmark_eigen();
    auto t2 = high_resolution_clock::now();
    std::cout << "Eigen sum = " << s1 << ", time = " 
              << duration_cast<milliseconds>(t2 - t1).count() << " ms\n";

    t1 = high_resolution_clock::now();
    float s2 = benchmark_vector_of_vector();
    t2 = high_resolution_clock::now();
    std::cout << "Vector of vector sum = " << s2 << ", time = " 
              << duration_cast<milliseconds>(t2 - t1).count() << " ms\n";

    t1 = high_resolution_clock::now();
    float s3 = benchmark_flattened_vector();
    t2 = high_resolution_clock::now();
    std::cout << "Flattened vector sum = " << s3 << ", time = " 
              << duration_cast<milliseconds>(t2 - t1).count() << " ms\n";

    return 0;
}
