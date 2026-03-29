#include <iostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <array>
#include <random>
#include <chrono>
#include <algorithm>

using namespace std;
using namespace std::chrono;

const vector<int> dims = {3, 9, 9, 9, 9, 9};     // 6D array
const int total_size = 1594323;
const int non_default_count = 7668;
const float default_value = 1000.0;
const int num_lookups = 1000;

// Block size (9x9x9)
const int BLOCK_LEN = 9;
const int BLOCK_SIZE = BLOCK_LEN * BLOCK_LEN * BLOCK_LEN;

int flatten_index(const vector<int>& idx, const vector<int>& dims) {
    int key = 0;
    for (size_t i = 0; i < dims.size(); ++i)
        key = key * dims[i] + idx[i];
    return key;
}

vector<int> generate_random_index(mt19937& rng) {
    vector<int> idx(dims.size());
    for (size_t i = 0; i < dims.size(); ++i) {
        uniform_int_distribution<int> dist(0, dims[i] - 1);
        idx[i] = dist(rng);
    }
    return idx;
}

// Blocked Sparse Array Access
int compute_block_id(const vector<int>& idx) {
    return idx[0] * (dims[1] * dims[2] * dims[3]) + idx[1] * dims[2] + idx[2];
}

int compute_offset(const vector<int>& idx) {
    return idx[3] * BLOCK_LEN * BLOCK_LEN + idx[4] * BLOCK_LEN + idx[5];
}

int main() {
    mt19937 rng(42);
    unordered_set<int> keys_set;

    vector<vector<int>> sparse_indices;

    while (sparse_indices.size() < non_default_count) {
        auto idx = generate_random_index(rng);
        int flat = flatten_index(idx, dims);
        if (keys_set.insert(flat).second) {
            sparse_indices.push_back(idx);
        }
    }

    vector<int> flat_keys;
    for (auto& idx : sparse_indices)
        flat_keys.push_back(flatten_index(idx, dims));

    // Random values
    uniform_real_distribution<float> val_dist(0.0, 999.0);
    vector<float> values(non_default_count);
    for (auto& v : values) v = val_dist(rng);

    // Dense array
    vector<float> dense_array(total_size, default_value);
    for (size_t i = 0; i < flat_keys.size(); ++i)
        dense_array[flat_keys[i]] = values[i];

    // Hash map
    unordered_map<int, float> hash_map;
    for (size_t i = 0; i < flat_keys.size(); ++i)
        hash_map[flat_keys[i]] = values[i];

    // Blocked sparse array
    unordered_map<int, array<float, BLOCK_SIZE>> block_map;
    for (size_t i = 0; i < sparse_indices.size(); ++i) {
        const auto& idx = sparse_indices[i];
        int block_id = compute_block_id(idx);
        int offset = compute_offset(idx);
        if (!block_map.count(block_id)) {
            block_map[block_id].fill(default_value);
        }
        block_map[block_id][offset] = values[i];
    }

    // Create 1000 hit and 1000 miss indices
    vector<int> hit_keys;
    vector<vector<int>> hit_indices;
    shuffle(sparse_indices.begin(), sparse_indices.end(), rng);
    for (int i = 0; i < num_lookups; ++i) {
        hit_indices.push_back(sparse_indices[i]);
        hit_keys.push_back(flatten_index(sparse_indices[i], dims));
    }

    vector<int> miss_keys;
    vector<vector<int>> miss_indices;
    while (miss_keys.size() < num_lookups) {
        auto idx = generate_random_index(rng);
        int flat = flatten_index(idx, dims);
        if (keys_set.count(flat) == 0) {
            miss_keys.push_back(flat);
            miss_indices.push_back(idx);
        }
    }

    float sum = 0;

    // Dense Array Hit
    auto start = high_resolution_clock::now();
    for (int key : hit_keys)
        sum += dense_array[key];
    auto end = high_resolution_clock::now();
    double dense_hit_time = duration<double>(end - start).count();

    // Dense Array Miss
    start = high_resolution_clock::now();
    for (int key : miss_keys)
        sum += dense_array[key];
    end = high_resolution_clock::now();
    double dense_miss_time = duration<double>(end - start).count();

    // Hash Map Hit
    start = high_resolution_clock::now();
    for (int key : hit_keys)
        sum += hash_map.count(key) ? hash_map[key] : default_value;
    end = high_resolution_clock::now();
    double map_hit_time = duration<double>(end - start).count();

    // Hash Map Miss
    start = high_resolution_clock::now();
    for (int key : miss_keys)
        sum += hash_map.count(key) ? hash_map[key] : default_value;
    end = high_resolution_clock::now();
    double map_miss_time = duration<double>(end - start).count();

    // Blocked Sparse Array Hit
    start = high_resolution_clock::now();
    for (const auto& idx : hit_indices) {
        int bid = compute_block_id(idx);
        int offset = compute_offset(idx);
        sum += block_map.count(bid) ? block_map[bid][offset] : default_value;
    }
    end = high_resolution_clock::now();
    double block_hit_time = duration<double>(end - start).count();

    // Blocked Sparse Array Miss
    start = high_resolution_clock::now();
    for (const auto& idx : miss_indices) {
        int bid = compute_block_id(idx);
        int offset = compute_offset(idx);
        sum += block_map.count(bid) ? block_map[bid][offset] : default_value;
    }
    end = high_resolution_clock::now();
    double block_miss_time = duration<double>(end - start).count();

    // Output
    cout << "Benchmark results (seconds for 1000 lookups):\n";
    cout << "-------------------------------------------------\n";
    cout << "Dense Array - Hit:   " << dense_hit_time << " s\n";
    cout << "Dense Array - Miss:  " << dense_miss_time << " s\n";
    cout << "Hash Map    - Hit:   " << map_hit_time << " s\n";
    cout << "Hash Map    - Miss:  " << map_miss_time << " s\n";
    cout << "Blocked Map - Hit:   " << block_hit_time << " s\n";
    cout << "Blocked Map - Miss:  " << block_miss_time << " s\n";

    return 0;
}