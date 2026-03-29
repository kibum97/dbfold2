#include "define.h"
#include "globals.h"
#include "misc_util.h"
Float rot_mat[3][3]; /* this is used in vector.h

                      */
#include "atom.h"
#include "getrms.h" /* getrms.h was moved up from includes at the bottom because it provides struct backbone */
#include "mpi.h"
#include "vector.h"

/* global dummy variables */
FILE *DATA, *DATA_OUT, *STATUS;

/* static data structures */
Float     *radii;

struct cell ***the_matrix;
long int       seed;


float          prob_ang[20][100];
float          deviation_ang[20][100][4];
int            no_chi_list[20];





// int max_N_cysteines=100; //Ideally would figure out the number of cysteines by getting the number
// of unique residues in the cons traints list..for now assume no more than 100

/* hydrogen bonding */
struct hydrogen_bond {
    float max_D2;
    float min_D2;
    float D_int;

    struct atom *donor;
    struct atom *donor2;
    struct atom *donor3;
    struct atom *donor4;
    struct atom *donor5;
    struct atom *donor6;
    struct atom *donor7;
    struct atom *acceptor;
    struct atom *acceptor2;
    struct atom *acceptor3;
    struct atom *acceptor4;
    struct atom *acceptor5;
    struct atom *acceptor6;
    struct atom *acceptor7;
    struct atom *acceptor8;
};
short                **type_contacts;
short             *****sct_E;
short                 *hbond_E;
struct hydrogen_bond **hbonds;
int                    n_hbond_int;
long int               NO_HBOND;
struct pair           *hbond_pair;
int                    total_hbond_pairs;
Float                  prev_E_hbond, E_hbond, dE_hbond;
float                  d_memory;

int   SEQ_DEP_HB;
float seq_hb[3][20][20];
int   helix_sheet;

char secstr[MAXSEQUENCE];
char path_dir[500]; /* 300>500 */

/* contact matrices for debugging */
#if DEBUG
unsigned char **debug_contacts;
unsigned char **debug_dcontacts;
short           debug_ncontacts;
unsigned char **debug_clashes;
short           debug_nclashes;
#endif



/* VZ:
 * The following two used to be local vars in fold.h, and struct_f2 also local to move.h
 * struct_f1 holds reference structure backbone.
 * I've moved them to here because I'd like the functions in move.h to access them
 * Yes even more global namespace pollution unfortunately.
 * A struct backbone holds five struct vector's named N, CA, C, CB, O
 * MAXSEQUENCE is macro for 350
 */
struct backbone struct_f1[MAXSEQUENCE], struct_f2[MAXSEQUENCE];
// struct backbone constraint_struct_f1[MAXSEQUENCE], constraint_struct_f2[MAXSEQUENCE];
Float Rg;                /* radius of gyration */
Float sc_rms, native_rms, new_rms, rms_Emin, rms_RMSDmin;
int   first_atom_res;
int  *helix;

/* potential globals */
Float **potential;
Float   mu, hydrogen_bond;

/* globals for move.h */
long int move_cycle = 0;
short           total_triplets;

short        total_pairs, total_pairs2;
short        TOTAL_SINGLE_LOOP_MOVES, TOTAL_DOUBLE_LOOP_MOVES, TOTAL_TRIPLE_LOOP_MOVES;
float        jacobi_before, jacobi_after;
int          soln_no_before, soln_no_after;
int          total_ntorsions;

/* alignment data structures */
struct align_cutoff {
    Float a;
    Float b;
}             **align_con_dist;
Float         **align_hard_core;

int             nalign, struct_natoms, struct_nresidues;
short           struct_ncontacts, struct_nclashes;
struct atom    *struct_native;
struct residue *struct_residue;
int            *seq_to_struct, *struct_to_seq;
int            *map_to_seq, *map_to_struct;
struct segment {
    int a;
    int b;
} str_segment[100], seq_segment[100];
int  nseg;
char amino_letters[20] = "ARNDCQEGHILKMFPSTWYV";
/* ----------------------- */

/*Bias*/
double  k_bias;

short **substructures;  // AB

int     umbrella;               /*Is umbrella simulation turned on?*/
int     number_of_contacts_max; /*Number of contacts setpoint for first node*/
int     contacts_step;          /*By how much does stepoint decrease for each node? */

int    new_natives;

int    n_substructures;
int   *substructure_sizes;          // AB

/* program option variables */
int           PRINT_PDB;
long int      MC_STEPS, MC_INIT_STEPS, MC_ANNEAL_STEPS;
int           SIDECHAIN_MOVES, MC_PRINT_STEPS, MC_PDB_PRINT_STEPS, MC_REPLICA_STEPS, MAX_EXCHANGE;
Float         MC_TEMP;
Float         MC_TEMP_MIN;
Float         TEMP_STEP;
int           NODES_PER_TEMP;
Float         STEP_SIZE;
Float         ALPHA, LAMBDA, beta;
Float         NON_NATIVE_REPULSION, NON_SPECIFIC_ENERGY, NATIVE_ATTRACTION;
Float         weight_clash, weight_potential, weight_rms;
Float         weight_hbond;
Float         LATTICE_SIZE;
Float         USE_GLOBAL_BB_MOVES;

Float         YANG_SCALE;
int           USE_SIDECHAINS, USE_ROTAMERS, USE_ROT_PROB, NO_NEW_CLASHES;
int           SKIP_LOCAL_CONTACT_RANGE, SKIP_BB_CONTACT_RANGE;
Float         SIDECHAIN_NOISE;
unsigned char MATRIX_SIZE, HALF_MATRIX_SIZE;
short         READ_POTENTIAL, READ_DENSITY, USE_GO_POTENTIAL;


Float  dE_tor, prev_E_tor, E_tor;
Float  dE_sct, prev_E_sct, E_sct;
Float  dE_aro, prev_E_aro, E_aro;
Float  dE_Rg, prev_E_Rg, E_Rg;
double a_PCA[MAXSEQUENCE], a_bCA[MAXSEQUENCE], phim[MAXSEQUENCE], psim[MAXSEQUENCE];
double cur_phi[MAXSEQUENCE], cur_psi[MAXSEQUENCE];
short  torsion_E[MAXSEQUENCE][6][6][6][6];
short  aromatic_E[9];
int    Naromatic;
int    Res_aromatic[MAXSEQUENCE];
float  cluster_phi[20][NOCLUSTERS];
float  cluster_psi[20][NOCLUSTERS];

/* parameter file names */

/* vzhao: I increased the length limits to 500 from 100/150 for these: */

char substructure_path[500];  // AB
/* file names added 12DEC02 IAH */
char min_etot_file[500], min_drms_file[500];


/* rotation data structures */
Float rot_mat_00, rot_mat_10, rot_mat_20, rot_mat_01, rot_mat_11, rot_mat_21, rot_mat_02,
    rot_mat_12, rot_mat_22;
short  **rotate_natoms;          /* number of atoms affected by the rotation */
short ***rotate_atom;            /* atoms affected by the rotation */
short  **loop_rotate_natoms;     /* number of atoms affected by the rotation */
short ***loop_rotate_atoms;      /* atoms affected by the rotation */
short  **loop_int_rotate_natoms; /* number of atoms affected by the rotation */
short ***loop_int_rotate_atoms;  /* atoms affected by the rotation */
char  ***loop_not_rotated;
char  ***not_rotated;
short  **rotate_sidechain_natoms;
short ***rotate_sidechain_atom;
char  ***sidechain_not_rotated;
short ***sidechain_torsion;
short   *all_rotated_atoms; /* pointer to atoms rotated by the current step */
short    all_rotated_natoms;
short    yang_rotated_natoms;
short   *yang_rotated_atoms;
char    *yang_not_rotated;

MPI_Comm   mpi_world_comm;
MPI_Status mpi_status;
int        nprocs, myrank, sel_num, ierr;
int  current_replica; /* Added 2021-07-12: to track a set of coordinates across replica space */
int *accepted_replica, *rejected_replica, *replica_index;
int *attempted_replica; /* VZ added 2-18 */
// int **accepted_replica, **rejected_replica, *replica_index; //AB changed to 2D array
float *Enode, *Tnode;
int   *Nnode, *Cnode;
float  delta_E, delta_T, delta_N, delta_all;
float *buf_in, *buf_out;

/*Things defined by AB related to knowledge-based moves*/



/*Added by AB to allow myrank to be offset, for instnace when running multi-core unfolding
 * simulations without MPI...by default, it is not*/
int my_rank_offset = 0;

/* includes */
#include "globals.h"

#include "align.h"
#include "constraint.h" /*AB on 11/26/19*/
#include "contacts.h"
#include "debug.h"
#include "energy.h"
#include "fold.h"
#include "hbonds.h"
#include "init.h"
#include "jac_local.h"
#include "lattice_util.h"
#include "loop.h"
#include "move.h"
#include "mu_potential.h"
#include "pdb_util.h"
#include "protein_util.h"
#include "rotate.h"
#include "tripep_closure.h"
#include "update.h"
