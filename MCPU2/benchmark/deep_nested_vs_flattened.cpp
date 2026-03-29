#include <iostream>
#include <chrono>

using namespace std;
using namespace std::chrono;

const int D1 = 6, D2 = 6, D3 = 6, D4 = 20, D5 = 20, D6 = 20;
const int N = D1 * D2 * D3 * D4 * D5 * D6;  // 864,000

template <typename F>
double benchmark(const string& label, F&& func) {
    auto start = high_resolution_clock::now();
    func();
    auto end = high_resolution_clock::now();
    double elapsed = duration<double>(end - start).count();
    cout << label << ": " << elapsed << " s" << endl;
    return elapsed;
}

int flatten(int i1, int i2, int i3, int i4, int i5, int i6) {
    return ((((i1 * D2 + i2) * D3 + i3) * D4 + i4) * D5 + i5) * D6 + i6;
}

int main() {
    cout << "Benchmarking deeply nested vs flattened array (6×6×6×20×20×20 = 864,000 elements)\n\n";

    // 1. Deeply nested array
    static float nested[D1][D2][D3][D4][D5][D6] = {};

    benchmark("Nested Fill", [&]() {
        for (int i1 = 0; i1 < D1; ++i1)
        for (int i2 = 0; i2 < D2; ++i2)
        for (int i3 = 0; i3 < D3; ++i3)
        for (int i4 = 0; i4 < D4; ++i4)
        for (int i5 = 0; i5 < D5; ++i5)
        for (int i6 = 0; i6 < D6; ++i6)
            nested[i1][i2][i3][i4][i5][i6] = static_cast<float>((i1 + i2 + i3 + i4 + i5 + i6) % 100);
    });

    float sum_nested = 0;
    benchmark("Nested Sum", [&]() {
        for (int i1 = 0; i1 < D1; ++i1)
        for (int i2 = 0; i2 < D2; ++i2)
        for (int i3 = 0; i3 < D3; ++i3)
        for (int i4 = 0; i4 < D4; ++i4)
        for (int i5 = 0; i5 < D5; ++i5)
        for (int i6 = 0; i6 < D6; ++i6)
            sum_nested += nested[i1][i2][i3][i4][i5][i6];
    });

    // 2. Flattened array
    float* flat = new float[N];

    benchmark("Flattened Fill", [&]() {
        for (int i1 = 0; i1 < D1; ++i1)
        for (int i2 = 0; i2 < D2; ++i2)
        for (int i3 = 0; i3 < D3; ++i3)
        for (int i4 = 0; i4 < D4; ++i4)
        for (int i5 = 0; i5 < D5; ++i5)
        for (int i6 = 0; i6 < D6; ++i6)
            flat[flatten(i1, i2, i3, i4, i5, i6)] = static_cast<float>((i1 + i2 + i3 + i4 + i5 + i6) % 100);
    });

    float sum_flat = 0;
    benchmark("Flattened Sum", [&]() {
        for (int i1 = 0; i1 < D1; ++i1)
        for (int i2 = 0; i2 < D2; ++i2)
        for (int i3 = 0; i3 < D3; ++i3)
        for (int i4 = 0; i4 < D4; ++i4)
        for (int i5 = 0; i5 < D5; ++i5)
        for (int i6 = 0; i6 < D6; ++i6)
            sum_flat += flat[flatten(i1, i2, i3, i4, i5, i6)];
    });

    cout << "\nCheck: sum_nested = " << sum_nested << ", sum_flat = " << sum_flat << "\n";

    delete[] flat;
    return 0;
}
