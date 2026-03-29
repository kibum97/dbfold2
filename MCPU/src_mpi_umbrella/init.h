#ifndef INIT_H
#define INIT_H

#include "atom.h"

extern FILE       *DATA;
extern char native_file[500], structure_file[500], native_directory[1000], triplet_file[500],
    sctorsion_file[500], sec_str_file[500], template_file[500], pdb_out_file[500],
    amino_data_file[500], rotamer_data_file[500], potential_file[500], atom_type_file[500],
    helicity_data[500], hydrogen_bonding_data[500];
extern char  substructure_path[500], alignment_file[500];
extern char std_file[500], std_prefix[500]; /* vzhao: added these here from backbone.c */
extern char aromatic_file[500];
extern char          PROTEIN_NAME[100];
extern int         PRINT_PDB, SIDECHAIN_MOVES, MC_PRINT_STEPS, MC_PDB_PRINT_STEPS, NODES_PER_TEMP;
extern int         SKIP_LOCAL_CONTACT_RANGE, SKIP_BB_CONTACT_RANGE, my_rank_offset, USE_SIDECHAINS;
extern int         NO_NEW_CLASHES, USE_ROTAMERS, READ_POTENTIAL, USE_GO_POTENTIAL, MC_REPLICA_STEPS;
extern int         MAX_EXCHANGE, umbrella, number_of_contacts_max, contacts_step, min_seq_sep;
extern long int    MC_STEPS, MC_ANNEAL_STEPS;
extern float       STEP_SIZE, SIDECHAIN_NOISE, MC_TEMP_MIN, TEMP_STEP, ALPHA, LAMBDA;
extern float       NATIVE_ATTRACTION, NON_NATIVE_REPULSION, USE_ROT_PROB, SEQ_DEP_HB;
extern float hydrogen_bond, weight_rms, rmsd_constraint, NON_SPECIFIC_ENERGY, USE_GLOBAL_BB_MOVES;
extern double YANG_MOVE, YANG_SCALE, contact_calpha_cutoff;
extern double weight_clash;

extern struct atom_type atom_type_list[200];
extern struct amino  *amino_acids;

void SetProgramOptions();
void InitializeProtein();
void DetermineTriplets();
void SetRadii();
void SetHardCore();
void SetContactDistance();
void ReadNative();
void GetResidueInfo(struct atom *, struct residue *, int, int);
void GetPhiPsi(struct atom *, struct residue *, int);
void GetChi();
void InitializeData();
void ReadSidechainTorsionData();
void InitializeSidechainRotationData();
void InitializeBackboneRotationData();
void CheckCorrelation(struct contact_data **, struct atom *, struct residue *, int);
void ReadPotential();
void ReadDistPotential();
int  SkipSelf(int, int, struct atom *, struct residue *);
int  SkipNeighbors(int, int, struct atom *, struct residue *);
void TurnOffNativeClashes(int);
void ReadAlignment();
void SetupAlignmentStructure();
void SetupAlignmentPotential();
void ReadAvgChis(void);
void ReadHelicityData();
void InitializeHydrogenBonding();
int  MatchAtomname(char *);

#endif
