#pragma once
#include <string_view>
#include <array>

struct HBondTopologyConfig {
    // string_view is much faster than string for compile-time constants in C++17/20
    std::string_view res1;
    std::string_view res2;
    
    // Hardcode the sizes to 7 and 8
    std::array<std::string_view, 7> donors;
    std::array<int, 7> donor_offsets;      
    
    std::array<std::string_view, 8> acceptors;
    std::array<int, 8> acceptor_offsets;   
    
    int skip;
    float min_D;
    float max_D;
    float D_int;
};

inline constexpr HBondTopologyConfig DEFAULT_HBOND_TOPOLOGY = {
    "XXX", "XXX",
    // Donors & Offsets (Exactly 7)
    {"N", "CA", "C", "C", "CA", "N", "CA"},
    {0, 0, -1, 0, -1, -1, 1},
    // Acceptors & Offsets (Exactly 8)
    {"O", "C", "CA", "N", "N", "CA", "C", "CA"},
    {0, 0, 0, 1, 0, 1, 1, -1},
    // Parameters
    4, 0.0f, 5.0f, 0.1f
};