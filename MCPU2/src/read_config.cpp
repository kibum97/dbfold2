#include "read_config.h"
#include "logger.h"
#include <yaml-cpp/yaml.h>
#include <iostream>
#include <fstream>
#include <string>

// Define external variables
FILE *DATA;
std::string native_file, structure_file, native_directory, template_file;
std::string substructure_path, alignment_file, triplet_file, sctorsion_file;
std::string sec_str_file, amino_data_file, atom_type_file, rotamer_data_file;
std::string pdb_out_file, std_prefix, potential_file, helicity_data;
std::string hydrogen_bonding_data, hbond_file, aromatic_file, PROTEIN_NAME;
int PRINT_PDB, SIDECHAIN_MOVES, MC_PRINT_STEPS, MC_PDB_PRINT_STEPS, NODES_PER_TEMP;
int SKIP_LOCAL_CONTACT_RANGE, SKIP_BB_CONTACT_RANGE, my_rank_offset, USE_SIDECHAINS;
int NO_NEW_CLASHES, USE_ROTAMERS, READ_POTENTIAL, USE_GO_POTENTIAL, MC_REPLICA_STEPS;
int MAX_EXCHANGE, umbrella, number_of_contacts_max, contacts_step, min_seq_sep;
int USE_CLUSTER, MAX_CLUSTERSTEP, CLUSTER_MOVE;
long int MC_STEPS, MC_ANNEAL_STEPS;
float STEP_SIZE, SIDECHAIN_NOISE, MC_TEMP_MIN, TEMP_STEP, ALPHA, LAMBDA;
float NATIVE_ATTRACTION, NON_NATIVE_REPULSION, USE_ROT_PROB, SEQ_DEP_HB, weight_clash;
float hydrogen_bond, weight_rms, rmsd_constraint, NON_SPECIFIC_ENERGY, USE_GLOBAL_BB_MOVES;
float YANG_MOVE, YANG_SCALE, contact_calpha_cutoff;

void SetProgramOptions(const std::string &cfg_file) {
    YAML::Node config = YAML::LoadFile(cfg_file);

    // Read native protein data
    native_file = config["native_protein_data"]["NATIVE_FILE"].as<std::string>();
    structure_file = config["native_protein_data"]["STRUCTURE_FILE"].as<std::string>();
    native_directory = config["native_protein_data"]["NATIVE_DIRECTORY"].as<std::string>();
    template_file = config["native_protein_data"]["TEMPLATE_FILE"].as<std::string>();
    alignment_file = config["native_protein_data"]["ALIGNMENT_FILE"].as<std::string>();
    pdb_out_file = config["native_protein_data"]["PDB_OUT_FILE"].as<std::string>();
    PROTEIN_NAME = config["native_protein_data"]["PROTEIN_NAME"].as<std::string>();

    // Read potential parameters
    NO_NEW_CLASHES = config["potential_parameters"]["NO_NEW_CLASHES"].as<int>();
    READ_POTENTIAL = config["potential_parameters"]["READ_POTENTIAL"].as<int>();
    USE_GO_POTENTIAL = config["potential_parameters"]["USE_GO_POTENTIAL"].as<int>();
    weight_clash = config["potential_parameters"]["CLASH_WEIGHT"].as<float>();
    weight_rms = config["potential_parameters"]["RMS_WEIGHT"].as<float>();
    hydrogen_bond = config["potential_parameters"]["HYDROGEN_BOND"].as<float>();
    NATIVE_ATTRACTION = config["potential_parameters"]["NATIVE_ATTRACTION"].as<float>();
    NON_NATIVE_REPULSION = config["potential_parameters"]["NON_NATIVE_REPULSION"].as<float>();
    NON_SPECIFIC_ENERGY = config["potential_parameters"]["NON_SPECIFIC_ENERGY"].as<float>();

    // Read contact definition
    SKIP_LOCAL_CONTACT_RANGE = config["contact_definition"]["SKIP_LOCAL_CONTACT_RANGE"].as<int>();
    SKIP_BB_CONTACT_RANGE = config["contact_definition"]["SKIP_BB_CONTACT_RANGE"].as<int>();

    // Read Monte-Carlo parameters
    rmsd_constraint = config["monte_carlo_parameters"]["CONSTRAINT_RMSD"].as<float>();

    // Read parameter files
    triplet_file = config["parameter_files"]["TRIPLET_ENERGY_FILE"].as<std::string>();
    sctorsion_file = config["parameter_files"]["SIDECHAIN_TORSION_FILE"].as<std::string>();
    sec_str_file = config["parameter_files"]["SECONDARY_STRUCTURE_FILE"].as<std::string>();
    amino_data_file = config["parameter_files"]["AMINO_DATA_FILE"].as<std::string>();
    rotamer_data_file = config["parameter_files"]["ROTAMER_DATA_FILE"].as<std::string>();
    atom_type_file = config["parameter_files"]["ATOM_TYPE_FILE"].as<std::string>();
    helicity_data = config["parameter_files"]["HELICITY_DATA"].as<std::string>();
    hydrogen_bonding_data = config["parameter_files"]["HYDROGEN_BONDING_DATA"].as<std::string>();
    hbond_file = config["parameter_files"]["HYDROGEN_BOND_FILE"].as<std::string>();
    potential_file = config["parameter_files"]["POTENTIAL_DATA"].as<std::string>();
    aromatic_file = config["parameter_files"]["AROMATIC_FILE"].as<std::string>();
}

void LogConfigurationDetails() {
    Logger::logf(INFO, "---GENERAL---");
    Logger::logf(INFO, "  VERSION: flat-bottom harmonic disulfide potential for all distances");
    if (CLUSTER_MOVE && USE_CLUSTER) {
        Logger::logf(INFO, "WARNING! KNOWLEDGE MOVES ARE IMPLEMENTED! THESE DO NOT SATISFY DETAILED BALANCE!  \n ");
    }
    Logger::logf(INFO, "Note that if any secondary structure character is anything except C, then the value of USE_CLUSTER specified in the cfg file will be overwritten by a default value (see move.h). \n This may lead to the implementation of (detailed-balance violating) knowledge moves!  \n");
    Logger::logf(INFO, "Starting configuration: %s", native_file.c_str());
    Logger::logf(INFO, "Target configuration:   %s", structure_file.c_str());
    Logger::logf(INFO, "Template Information:   %s", template_file.c_str());
    Logger::logf(INFO, "Potential data:         %s", potential_file.c_str());
    Logger::logf(INFO, "Atom typing scheme:     %s", atom_type_file.c_str());
    Logger::logf(INFO, "Read potential:         %d", READ_POTENTIAL);
    /*
    Logger::logf(INFO, "Skip local contact:     %d\n", SKIP_LOCAL_CONTACT_RANGE);
    Logger::logf(INFO, "Skip bb contact:        %d\n\n", SKIP_BB_CONTACT_RANGE);
    Logger::logf(INFO, "Alpha:               %7.2f\n", ALPHA);
    Logger::logf(INFO, "Beta:                %7.2f\n", beta);
    Logger::logf(INFO, "Lambda:              %7.2f\n\n", LAMBDA);
    Logger::logf(INFO, "STEP_SIZE:           %7.2f\n\n", STEP_SIZE * rad2deg);
    Logger::logf(INFO, "YANG_MOVE:           %7.2f\n", YANG_MOVE);
    Logger::logf(INFO, "YANG_SCALE:          %7.2f\n\n", YANG_SCALE);
    Logger::logf(INFO, "CLUSTER_MOVE:        %7.2f\n", CLUSTER_MOVE);
    Logger::logf(INFO, "USE_CLUSTER:         %7.2f\n", USE_CLUSTER);
    Logger::logf(INFO, "NOCLUSTERS:          %4d\n", NOCLUSTERS);
    Logger::logf(INFO, "CLUSTER_NOISE:       %7.2f\n\n", CLUSTER_NOISE);
    Logger::logf(INFO, "beta_favor:          %7.2f\n", beta_favor);
    Logger::logf(INFO, "HB_CUTOFF:           %7.2f\n", HB_CUTOFF);
    Logger::logf(INFO, "HB_INNER:            %7.2f\n", HB_INNER);
    Logger::logf(INFO, "HB_PENALTY:          %7.2f\n", HB_PENALTY);
    Logger::logf(INFO, "SEQ_DEP_HB:          %7d\n", SEQ_DEP_HB);
    Logger::logf(INFO, "RDTHREE_CON:         %7.2f\n", RDTHREE_CON);
    Logger::logf(INFO, "AROMATIC_DISTANCE:   %7.2f\n\n", AROMATIC_DISTANCE);
    Logger::logf(INFO, "USE_ROTAMERS         %4d\n", USE_ROTAMERS);
    Logger::logf(INFO, "USE_ROT_PROB         %4d\n\n", USE_ROT_PROB);
    Logger::logf(INFO, "SIDECHAIN_MOVES      %4d\n\n", SIDECHAIN_MOVES);
    Logger::logf(INFO, "MC_REPLICA_STEPS:        %7d\n", MC_REPLICA_STEPS);
    Logger::logf(INFO, "exchanges in each step:  %7d\n\n", MAX_EXCHANGE);
    Logger::logf(INFO, "MC_STEPS:            %10ld\n", MC_STEPS);
    Logger::logf(INFO, "MC_ANNEAL_STEPS:     %10ld\n", MC_ANNEAL_STEPS);
    Logger::logf(INFO, "MC_PRINT_STEPS:      %10d\n", MC_PRINT_STEPS);
    Logger::logf(INFO, "MC_PDB_PRINT_STEPS:  %10d\n\n", MC_PDB_PRINT_STEPS);
    */
}