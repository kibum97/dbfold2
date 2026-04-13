#include "ConfigLoader.hpp"
#include "globals.h" 

void ConfigLoader::read(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open config file: " + filename);
    }
    config = nlohmann::json::parse(file);
}

void ConfigLoader::loadSystem(System& sys) {
    // Instead of auto& pp = config["..."], we use the helpers directly
    sys.USE_GO_POTENTIAL      = get<int>("potential_parameters", "USE_GO_POTENTIAL");
    sys.weight_clash          = get<float>("potential_parameters", "CLASH_WEIGHT");
    sys.weight_rms            = get<float>("potential_parameters", "RMS_WEIGHT");
    sys.NATIVE_ATTRACTION     = get<float>("potential_parameters", "NATIVE_ATTRACTION");
    sys.NON_NATIVE_REPULSION  = get<float>("potential_parameters", "NON_NATIVE_REPULSION");
    sys.NON_SPECIFIC_ENERGY   = get<float>("potential_parameters", "NON_SPECIFIC_ENERGY");
}

void ConfigLoader::loadSimulation(Simulation& sim) {
    // Required values use 'get'
    sim.native_file      = get<std::string>("native_protein_data", "NATIVE_FILE");
    sim.structure_file   = get<std::string>("native_protein_data", "STRUCTURE_FILE");
    sim.PROTEIN_NAME     = get<std::string>("native_protein_data", "PROTEIN_NAME");
    
    // Optional values use 'get_or'
    sim.native_directory = get_or<std::string>("native_protein_data", "NATIVE_DIRECTORY", "None");
    
    sim.SKIP_LOCAL_CONTACT_RANGE = get<int>("contact_definition", "SKIP_LOCAL_CONTACT_RANGE");
}

void ConfigLoader::loadIntegrator(MCIntegrator& integrator) {
    integrator.YANG_MOVE      = get<float>("monte_carlo_parameters", "YANG_MOVE");
    integrator.YANG_SCALE     = get<int>("monte_carlo_parameters", "YANG_SCALE");
    integrator.USE_SIDECHAINS = get<int>("monte_carlo_parameters", "USE_SIDECHAINS");
}