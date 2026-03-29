#include <iostream>
#include <vector>
#include <array>
#include <chrono>
#include <numeric>

using namespace std;
using namespace std::chrono;

const size_t N = 10'000'000;

template <typename T>
double benchmark_fill_and_sum(T& container, const string& name) {
    // Fill
    auto start = high_resolution_clock::now();
    for (size_t i = 0; i < N; ++i)
        container[i] = static_cast<float>(i % 100);
    // Sum
    float sum = 0;
    for (size_t i = 0; i < N; ++i)
        sum += container[i];
    auto end = high_resolution_clock::now();

    duration<double> elapsed = end - start;
    cout << name << " → Sum: " << sum << ", Time: " << elapsed.count() << " s" << endl;
    return elapsed.count();
}

int main() {
    cout << "Benchmarking fill + sum for " << N << " floats...\n\n";

    // 1. C-style array (on heap to avoid stack overflow)
    float* c_array = new float[N];
    benchmark_fill_and_sum(c_array, "C-style array");
    delete[] c_array;

    // 2. std::array (on heap using pointer workaround)
    auto heap_array = new std::array<float, N>;
    benchmark_fill_and_sum(*heap_array, "std::array (heap)");
    delete heap_array;

    // 3. std::vector
    vector<float> vec(N);
    benchmark_fill_and_sum(vec, "std::vector");

    return 0;
}
