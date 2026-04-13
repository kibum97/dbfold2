#include "globals.h"

struct Topology top;
struct System sys;
struct MCIntegrator integrator;
struct Context ctx;
struct Simulation sim;

/* Global Variable Definitions */
int   natoms = 0, nresidues = 0, nprocs = 0, myrank = 0, sel_num = 0, ierr = 0, current_replica = 0;
FILE *STATUS = NULL;
MPI_Comm mpi_world_comm;

struct atom *native = NULL, *prev_native = NULL, *orig_native = NULL, *native_Emin = NULL,
            *native_RMSDmin   = NULL;
struct contact_data    **data = NULL, **struct_data = NULL;
struct residue          *native_residue = NULL;
struct monte_carlo       mc             = {0};
Float                    rot_mat[3][3]  = {0};
short                    nclashes = 0, ncontacts = 0, delta_nclashes = 0, delta_contacts = 0;
struct cell           ***the_matrix = NULL;
struct pair             *ab = NULL, *cd = NULL;
int                      MAX_TYPES        = 0;
short                  **type_contacts    = NULL;
long int               **hard_core        = NULL;
struct cutoff          **contact_distance = NULL;
struct monte_carlo_flags mc_flags         = {0};
float                    USE_CLUSTER      = 0.0f;
int                      natom_type_list = 0, bb_O_type = 0, bb_N_type = 0, bb_OXT_type = 0;
struct atom             *temp_atom = NULL, *temp_prev_atom = NULL;
struct int_vector        temp_xyz_int = {0};
struct cell             *temp_cell = NULL, *temp_cell3 = NULL;
struct cell            **temp_cell_array = NULL;
unsigned char           *is_rotated      = NULL;
int      total_ntorsions = 0, total_pairs = 0, total_pairs2 = 0, total_hbond_pairs = 0;
char     secstr[MAXSEQUENCE] = {0};
long int distance            = 0;

/* Move related variables */
int TOTAL_TRIPLE_LOOP_MOVES = 0, TOTAL_SINGLE_LOOP_MOVES = 0, TOTAL_DOUBLE_LOOP_MOVES = 0;
struct triplet *residue_triplets     = NULL;
int             is_template[MAX_RES] = {0};
int             naccepted = 0, n_sidechain_accepted = 0;
int             nrejected = 0, nothers = 0, nomove = 0;
double          deg2rad            = 0.0;
int             move_cycle         = 0;
char         ***not_rotated        = NULL;
short        ***rotate_atom        = NULL;
short         **rotate_natoms      = NULL;
double        **cluster_phi        = NULL;
double        **cluster_psi        = NULL;
double         *cur_phi            = NULL;
double         *cur_psi            = NULL;
short        ***loop_rotate_atoms  = NULL;
short         **loop_rotate_natoms = NULL;
char         ***loop_not_rotated   = NULL;
int             no_chi_list[20]    = {0};
double        **prob_ang           = NULL;
int            *cur_rotamers = NULL, old_rotamer = 0;
short        ***sidechain_torsion         = NULL;
float           deviation_ang[20][100][4] = {0};
struct angles  *rotamer_angles            = NULL;
int             sidemovedone              = 0;
short          *all_rotated_atoms         = NULL;
short           all_rotated_natoms        = 0;
short         **rotate_sidechain_natoms   = NULL;
short        ***rotate_sidechain_atom     = NULL;
char         ***sidechain_not_rotated     = NULL;
int            *disulfide_pairs           = NULL;
int            *disulfide_pairs_attempt   = NULL;
int             soln_no_before = 0, jacobi_after = 0, jacobi_before = 0;
double          k_bias                          = 0.0;
int             diff_number_of_contacts_new     = 0;
int             diff_number_of_contacts_current = 0;
short          *orig_contactstring              = NULL;
int             new_natives                     = 0;
struct backbone struct_f1[MAXSEQUENCE] = {0}, struct_f2[MAXSEQUENCE] = {0};
int             sidechain_step              = 0;
int             backbone_accepted           = 0;
int             natives                     = 0;
int             number_of_contacts_setpoint = 0;
float           CLUSTER_MOVE                = 0.0f;
Float           YANG_MOVE                   = 0.0;

char constraint_file[500]      = {0};
char rmsd_constraint_file[500] = {0};
char movable_region_file[500]  = {0};
int  nalign                    = 0;

struct pair *hbond_pair      = NULL;
Float        native_E        = 0.0;
Float       *radii           = NULL;
int          SEQ_DEP_HB      = 0;
int          MAX_CLUSTERSTEP = 0;
int         *replica_index = NULL, *accepted_replica = NULL, *rejected_replica = NULL,
    *attempted_replica  = NULL;
Float      LATTICE_SIZE = 0.0;
int        MATRIX_SIZE = 0, HALF_MATRIX_SIZE = 0;
short *****sct_E          = NULL;
int        total_triplets = 0;

/* Debugging Matrices */
unsigned char **debug_contacts  = NULL;
unsigned char **debug_dcontacts = NULL;
unsigned char **debug_clashes   = NULL;

/* Backbone Rotation Tracking */
short  **loop_int_rotate_natoms = NULL;
short ***loop_int_rotate_atoms  = NULL;

/* Energy and Simulation Variables */
Float E = 0.0, E_pot = 0.0, dE = 0.0, dE_pot = 0.0, prev_E = 0.0, prev_E_pot = 0.0;
Float E_tor = 0.0, E_sct = 0.0, E_aro = 0.0, E_Rg = 0.0, dE_tor = 0.0, dE_sct = 0.0, dE_aro = 0.0,
      dE_Rg      = 0.0;
Float prev_E_aro = 0.0, prev_E_sct = 0.0, prev_E_tor = 0.0;
float Emin = 0.0f, Emin_pot = 0.0f, Emin_hbond = 0.0f, Emin_tor = 0.0f, Emin_sct = 0.0f,
      Emin_aro = 0.0f, Emin_constraint = 0.0f;
Float    prev_E_hbond = 0.0, E_hbond = 0.0, dE_hbond = 0.0;
long int mcstep = 0, seed = 0;
Float  **potential = NULL;
Float    mu = 0.0, hydrogen_bond = 0.0;
double   new_rms = 0.0;

/* Hydrogen Bond Variables */
float                  seq_hb[3][20][20] = {0};
struct hydrogen_bond **hbonds            = NULL;
short                 *hbond_E           = NULL;
int                    helix_sheet       = 0;
long int               NO_HBOND          = 0;
float                  d_memory          = 0.0f;
short                  M = 0, N = 0, O = 0, P = 0;
short                 *A = NULL, *B = NULL;

/* Alignment Variables */
struct align_cutoff **align_con_dist  = NULL;
Float               **align_hard_core = NULL;
struct atom          *struct_native   = NULL;
struct residue       *struct_residue  = NULL;
int                   struct_natoms = 0, struct_nresidues = 0;
short                 struct_ncontacts = 0, struct_nclashes = 0;
int                  *seq_to_struct = NULL, *struct_to_seq = NULL;
int                  *map_to_seq = NULL, *map_to_struct = NULL;
struct segment        str_segment[100] = {0}, seq_segment[100] = {0};
int                   nseg  = 0;
int                  *helix = NULL;

/* Loop Variables */
int    n_movable_residues           = 0;
int    movable_residue_map[MAX_RES] = {0};
int    res_atomno[MAX_RES]          = {0};
short  yang_rotated_natoms          = 0;
short *yang_rotated_atoms           = NULL;
char  *yang_not_rotated             = NULL;
double a_PCA[MAXSEQUENCE] = {0}, a_bCA[MAXSEQUENCE] = {0}, phim[MAXSEQUENCE] = {0},
       psim[MAXSEQUENCE]                 = {0};
short torsion_E[MAXSEQUENCE][6][6][6][6] = {0};
short aromatic_E[9]                      = {0};
int   Naromatic                          = 0;
int   Res_aromatic[MAXSEQUENCE]          = {0};

/* Rotation variables */
Float rot_mat_00 = 0.0, rot_mat_10 = 0.0, rot_mat_20 = 0.0, rot_mat_01 = 0.0, rot_mat_11 = 0.0,
      rot_mat_21 = 0.0, rot_mat_02 = 0.0, rot_mat_12 = 0.0, rot_mat_22 = 0.0;

/* Debugging Variables */
#if DEBUG
short debug_ncontacts = 0;
short debug_nclashes  = 0;
#endif

// pdb utils variable
int first_atom_res = 0;

// fold variables
Float      native_rms = 0.0, rms_RMSDmin = 0.0, rms_Emin = 0.0;
long int   mcstep_Emin = 0, mcstep_RMSDmin = 0;
float      delta_E = 0.0f, delta_T = 0.0f, delta_N = 0.0f, delta_all = 0.0f;
MPI_Status mpi_status;
float      E_RMSDmin     = 0.0f;
float      E_RMSDmin_pot = 0.0f, E_RMSDmin_hbond = 0.0f, E_RMSDmin_tor = 0.0f, E_RMSDmin_sct = 0.0f,
      E_RMSDmin_aro = 0.0f, E_RMSDmin_constraint = 0.0f;

/* backbone variables */
char path_dir[500] = {0};

/* Constraint Variables */
int    Max_N_constraints = 0, N_constraints = 0;
int  **constraint_array = NULL;
float *distances = NULL, *constraint_weights = NULL, k_constraint = 0.0f, E_constraint = 0.0f,
      prev_E_constraint = 0.0f, dE_constraint = 0.0f, mean_constraint_distance = 0.0f;

/* MPI Buffers */
float *Enode = NULL, *Tnode = NULL, *buf_in = NULL, *buf_out = NULL;
int   *Nnode = NULL, *Cnode = NULL;

/* File Pointers */
FILE *DATA = NULL;

/* Strings and File Paths */
char native_file[500]           = {0};
char structure_file[500]        = {0};
char native_directory[1000]     = {0};
char triplet_file[500]          = {0};
char sctorsion_file[500]        = {0};
char sec_str_file[500]          = {0};
char template_file[500]         = {0};
char pdb_out_file[500]          = {0};
char amino_data_file[500]       = {0};
char rotamer_data_file[500]     = {0};
char potential_file[500]        = {0};
char atom_type_file[500]        = {0};
char helicity_data[500]         = {0};
char hydrogen_bonding_data[500] = {0};
char substructure_path[500]     = {0};
char alignment_file[500]        = {0};
char std_file[500]              = {0};
char std_prefix[500]            = {0};
char aromatic_file[500]         = {0};
char PROTEIN_NAME[100]          = {0};

/* Integer Flags and Counters */
int PRINT_PDB                = 0;
int SIDECHAIN_MOVES          = 0;
int MC_PRINT_STEPS           = 0;
int MC_PDB_PRINT_STEPS       = 0;
int NODES_PER_TEMP           = 0;
int SKIP_LOCAL_CONTACT_RANGE = 0;
int SKIP_BB_CONTACT_RANGE    = 0;
int my_rank_offset           = 0;
int USE_SIDECHAINS           = 0;
int NO_NEW_CLASHES           = 0;
int USE_ROTAMERS             = 0;
int MAX_EXCHANGE             = 0;
int umbrella                 = 0;
int number_of_contacts_max   = 0;
int contacts_step            = 0;
int min_seq_sep              = 0;
int USE_ROT_PROB             = 0;
int MC_REPLICA_STEPS         = 0;

/* Long Integers */
long int MC_STEPS        = 0;
long int MC_ANNEAL_STEPS = 0;
long int MC_INIT_STEPS   = 0;

/* Short Integers */
short READ_POTENTIAL   = 0;
short READ_DENSITY     = 0;
short USE_GO_POTENTIAL = 0;

/* Standard Floats and Doubles */
float rmsd_constraint = 0.0f;

/* Custom Typedef 'Float' Variables */
Float STEP_SIZE            = 0.0;
Float SIDECHAIN_NOISE      = 0.0;
Float MC_TEMP_MIN          = 0.0;
Float TEMP_STEP            = 0.0;
Float ALPHA                = 0.0;
Float LAMBDA               = 0.0;
Float MC_TEMP              = 0.0;
Float beta                 = 0.0;
Float NON_NATIVE_REPULSION = 0.0;
Float NON_SPECIFIC_ENERGY  = 0.0;
Float NATIVE_ATTRACTION    = 0.0;
Float weight_clash         = 0.0;
Float weight_potential     = 0.0;
Float weight_rms           = 0.0;
Float weight_hbond         = 0.0;
Float USE_GLOBAL_BB_MOVES  = 0.0;
Float YANG_SCALE           = 0.0;

/* Structs and Arrays */
struct atom_type atom_type_list[200] = {0};
struct amino    *amino_acids         = NULL;

FILE *fpdb = NULL;