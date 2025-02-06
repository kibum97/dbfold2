#ifndef READ_CONFIG_H
#define READ_CONFIG_H

#include <cstdio>
#include <string>

// Declare external variables
extern FILE *DATA;
extern std::string parameter_file;
extern std::string native_file, structure_file, native_directory, template_file;
extern std::string substructure_path, alignment_file, triplet_file, sctorsion_file;
extern std::string sec_str_file, amino_data_file, atom_type_file, rotamer_data_file;
extern std::string pdb_out_file, std_prefix, helicity_data;
extern std::string hydrogen_bonding_data, hbond_file, seq_hbond_file, aromatic_file, PROTEIN_NAME;
extern std::string all_triplet_file, all_sctorsion_file;
extern int PRINT_PDB, SIDECHAIN_MOVES, MC_PRINT_STEPS, MC_PDB_PRINT_STEPS, NODES_PER_TEMP;
extern int SKIP_LOCAL_CONTACT_RANGE, SKIP_BB_CONTACT_RANGE, my_rank_offset, USE_SIDECHAINS;
extern int NO_NEW_CLASHES, USE_ROTAMERS, READ_POTENTIAL, USE_GO_POTENTIAL, MC_REPLICA_STEPS;
extern int MAX_EXCHANGE, umbrella, number_of_contacts_max, contacts_step, min_seq_sep;
extern int USE_CLUSTER, MAX_CLUSTERSTEP, CLUSTER_MOVE;
extern long int MC_STEPS, MC_ANNEAL_STEPS;
extern float STEP_SIZE, SIDECHAIN_NOISE, MC_TEMP_MIN, TEMP_STEP, ALPHA, LAMBDA;
extern float NATIVE_ATTRACTION, NON_NATIVE_REPULSION, USE_ROT_PROB, SEQ_DEP_HB, weight_clash;
extern float hydrogen_bond, weight_rms, rmsd_constraint, NON_SPECIFIC_ENERGY, USE_GLOBAL_BB_MOVES;
extern float YANG_MOVE, YANG_SCALE, contact_calpha_cutoff;

// Function declaration
void SetProgramOptions(const std::string &cfg_file);
void LogConfigurationDetails();

#endif // READ_CONFIG_H