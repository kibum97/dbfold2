#pragma once

#include <nlohmann/json.hpp>
#include <fstream>
#include <string>
#include <stdexcept>

// Forward declarations of your structs to avoid including their headers here
struct Simulation;
struct System;
struct MCIntegrator;

class ConfigLoader {
public:
    nlohmann::json config;

    void read(const std::string& filename);

    // Templates MUST stay in the header
    template<typename T>
    T get(const std::string& cat, const std::string& key) const {
        return config.at(cat).at(key).get<T>();
    }

    template<typename T>
    T get_or(const std::string& cat, const std::string& key, T fallback) const {
        if (config.contains(cat) && config[cat].contains(key)) {
            return config[cat][key].get<T>();
        }
        return fallback;
    }

    // High-level mapping functions (defined in .cpp)
    void loadSystem(System& sys);
    void loadSimulation(Simulation& sim);
    void loadIntegrator(MCIntegrator& integrator);
};