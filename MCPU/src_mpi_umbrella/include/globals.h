#ifndef GLOBALS_H
#define GLOBALS_H

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "define.h"

/* Math constants */
#define pi 3.141592653589793238462643383279502884197e0
#define rad2deg 180.0e0 / pi

/* Struct definitions */
struct int_vector {
    long int x, y, z;
};
struct vector {
    float x, y, z;
};
struct pair {
    short a, b;
};

struct atom {
    struct vector     xyz;
    struct int_vector xyz_int;
    char              res[5], atomname[5];
    short         res_num, atomtype, smogtype, is_core, is_designed, is_sidechain, sec_structure;
    unsigned char X, Y, Z;
    struct cell  *matrix;
};

struct atom_type {
    char res_name[4];
    char atom_name[4];
    int  type_num;
};

struct amino {
    char  name[4];
    char  symbol[2];
    char  torsion[4][4][4];
    Float avg_angle[4][4][4][4][4];
    int   ntorsions;
    int   nrotamers;
    int   rotate_natoms[4];
    char  rotate_atom[4][10][4];
};

struct residue {
    char    res[4];
    int     amino_num, ntorsions, nrotamers;
    short   is_core, is_designed;
    Float   psi, phi, chi[4], native_chi[4], tmpchi[4];
    Float   avg_angle[4][4][4][4][4];
    short **rot_position;
    int     CA, N, C, O, CB, CG, CE1, CE2, CZ2, CZ3, atomnumber[40];
};

struct cell {
    char         X, Y, Z;
    struct cell *neighbors[27];
    short atom_list[MAX_CELL_ATOMS];  // AB: Each cell has this array called atom_list which has
                                      // MAX_CELL_ATOMS = 100 elements
    unsigned char natoms;
};

struct backbone {
    struct vector N;
    struct vector CA;
    struct vector C;
    struct vector CB;
    struct vector O;
};

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

struct monte_carlo {
    int backbone_move, sel_res_num, sel_rotamer, is_phi, torsion, loop_size, sel_triplet,
        selected[3];
    Float delta_angle[4], delta_phi_angle[3], delta_psi_angle[3];
};

struct monte_carlo_flags {
    unsigned char clashed : 1;
    unsigned char init : 1;
};

struct contact_data {
    unsigned char contacts : 1, delta_contacts : 1, clashes : 1, delta_clashes : 1, disulfide : 1;
    unsigned char check_contacts : 1, check_clashes : 1, check_hbond : 1;
    long int      hbond, delta_hbond;
    int           closehb, delta_closehb;
};

struct cutoff {
    long int a;
    long int b;
};

// Geometry related structs
struct triplet {
    short a;
    short b;
    short c;
};

struct angles {
    Float chis[100][4];
};

/* Global Variable Declarations (Extern Only) */
extern int                      natoms, nresidues, nprocs, myrank, sel_num, ierr, current_replica;
extern FILE                    *STATUS;
extern MPI_Comm                 mpi_world_comm;
extern struct atom             *native, *prev_native, *orig_native, *native_Emin, *native_RMSDmin;
extern struct contact_data    **data, **struct_data;
extern struct residue          *native_residue;
extern struct monte_carlo       mc;
extern Float                    rot_mat[3][3];
extern short                    nclashes, ncontacts, delta_nclashes, delta_contacts;
extern struct cell           ***the_matrix;
extern struct pair             *ab, *cd;
extern int                      MAX_TYPES;
extern short                  **type_contacts;
extern long int               **hard_core;
extern struct cutoff          **contact_distance;
extern struct monte_carlo_flags mc_flags;
extern float                    USE_CLUSTER;
extern int                      MAX_TYPES, natom_type_list, bb_O_type, bb_N_type, bb_OXT_type;
extern struct atom             *temp_atom, *temp_prev_atom;
extern struct int_vector        temp_xyz_int;
extern struct cell             *temp_cell, *temp_cell3;
extern struct cell            **temp_cell_array;
extern unsigned char           *is_rotated;
extern int                      min_seq_sep;
extern int                      total_ntorsions, total_pairs, total_pairs2, total_hbond_pairs;
extern char                     secstr[MAXSEQUENCE];
extern long int                 distance;

/* Move related variables */
extern int             TOTAL_TRIPLE_LOOP_MOVES, TOTAL_SINGLE_LOOP_MOVES, TOTAL_DOUBLE_LOOP_MOVES;
extern struct triplet *residue_triplets;
extern int             is_template[MAX_RES];
extern int             naccepted, n_sidechain_accepted;
extern int             nrejected, nothers, nomove;
extern double          deg2rad;
extern int             move_cycle;
extern char         ***not_rotated;
extern short        ***rotate_atom;
extern short         **rotate_natoms;
extern double        **cluster_phi;
extern double        **cluster_psi;
extern double         *cur_phi;
extern double         *cur_psi;
extern short        ***loop_rotate_atoms;
extern short         **loop_rotate_natoms;
extern char         ***loop_not_rotated;
extern int             no_chi_list[20];
extern double        **prob_ang;
extern int            *cur_rotamers, old_rotamer;
extern short        ***sidechain_torsion;
extern float           deviation_ang[20][100][4];
extern struct angles  *rotamer_angles;
extern int             sidemovedone;
extern short          *all_rotated_atoms;
extern short           all_rotated_natoms;
extern short         **rotate_sidechain_natoms;
extern short        ***rotate_sidechain_atom;
extern char         ***sidechain_not_rotated;
extern int            *disulfide_pairs;
extern int            *disulfide_pairs_attempt;
extern int             soln_no_before, jacobi_after, jacobi_before;
extern double          k_bias;
extern int             diff_number_of_contacts_new;
extern int             diff_number_of_contacts_current;
extern short          *orig_contactstring;
extern int             new_natives;
extern struct backbone struct_f1[MAXSEQUENCE], struct_f2[MAXSEQUENCE];
extern int             sidechain_step;
extern int             backbone_accepted;
extern int             natives;
extern int             number_of_contacts_setpoint;
extern float           CLUSTER_MOVE;
extern Float           YANG_MOVE;

extern char constraint_file[500];
extern char rmsd_constraint_file[500];
extern char movable_region_file[500];
extern int  nalign;

extern struct pair *hbond_pair;
extern Float        native_E;
extern Float       *radii;
extern int          SEQ_DEP_HB;
extern int          MAX_CLUSTERSTEP;
extern int         *replica_index, *accepted_replica, *rejected_replica, *attempted_replica;
extern Float        LATTICE_SIZE;
extern int          MATRIX_SIZE, HALF_MATRIX_SIZE;
extern short   *****sct_E;
extern int          total_triplets;

/* Debugging Matrices */
extern unsigned char **debug_contacts;
extern unsigned char **debug_dcontacts;
extern unsigned char **debug_clashes;

/* Backbone Rotation Tracking */
extern short  **loop_int_rotate_natoms;
extern short ***loop_int_rotate_atoms;

/* Energy and Simulation Variables */
extern Float    E, E_pot, dE, dE_pot, prev_E, prev_E_pot;
extern Float    E_tor, E_sct, E_aro, E_Rg, dE_tor, dE_sct, dE_aro, dE_Rg;
extern Float    prev_E_aro, prev_E_sct, prev_E_tor;
extern float    Emin, Emin_pot, Emin_hbond, Emin_tor, Emin_sct, Emin_aro, Emin_constraint;
extern Float    prev_E_hbond, E_hbond, dE_hbond;
extern long int mcstep, seed;
extern Float  **potential;
extern Float    mu, hydrogen_bond;
extern double   new_rms;

/* Hydrogen Bond Variables */
extern float                  seq_hb[3][20][20];
extern struct hydrogen_bond **hbonds;
extern short                 *hbond_E;
extern int                    helix_sheet;
extern long int               NO_HBOND;
extern float                  d_memory;
extern short                  M, N, O, P;
extern short                 *A, *B;

/* Initialization Variables */
extern FILE *DATA;
extern char  native_file[500], structure_file[500], native_directory[1000], triplet_file[500],
    sctorsion_file[500], sec_str_file[500], template_file[500], pdb_out_file[500],
    amino_data_file[500], rotamer_data_file[500], potential_file[500], atom_type_file[500],
    helicity_data[500], hydrogen_bonding_data[500];
extern char     substructure_path[500], alignment_file[500];
extern char     std_file[500], std_prefix[500]; /* vzhao: added these here from backbone.c */
extern char     aromatic_file[500];
extern char     PROTEIN_NAME[100];
extern int      PRINT_PDB, SIDECHAIN_MOVES, MC_PRINT_STEPS, MC_PDB_PRINT_STEPS, NODES_PER_TEMP;
extern int      SKIP_LOCAL_CONTACT_RANGE, SKIP_BB_CONTACT_RANGE, my_rank_offset, USE_SIDECHAINS;
extern int      NO_NEW_CLASHES, USE_ROTAMERS;
extern int      MAX_EXCHANGE, umbrella, number_of_contacts_max, contacts_step, min_seq_sep;
extern long int MC_STEPS, MC_ANNEAL_STEPS;
extern Float    STEP_SIZE, SIDECHAIN_NOISE, MC_TEMP_MIN, TEMP_STEP, ALPHA, LAMBDA;
extern int      USE_ROT_PROB;
extern float    hydrogen_bond, rmsd_constraint;
extern double   contact_calpha_cutoff;

extern struct atom_type atom_type_list[200];
extern struct amino    *amino_acids;

/* program option variables */
extern int      PRINT_PDB;
extern long int MC_STEPS, MC_INIT_STEPS, MC_ANNEAL_STEPS;
extern int      SIDECHAIN_MOVES, MC_PRINT_STEPS, MC_PDB_PRINT_STEPS, MC_REPLICA_STEPS, MAX_EXCHANGE;
extern Float    MC_TEMP;
extern Float    MC_TEMP_MIN;
extern Float    TEMP_STEP;
extern int      NODES_PER_TEMP;
extern Float    STEP_SIZE;
extern Float    ALPHA, LAMBDA, beta;
extern Float    NON_NATIVE_REPULSION, NON_SPECIFIC_ENERGY, NATIVE_ATTRACTION;
extern Float    weight_clash, weight_potential, weight_rms;
extern Float    weight_hbond;
extern Float    USE_GLOBAL_BB_MOVES;

extern Float YANG_SCALE;
extern int   USE_SIDECHAINS, USE_ROTAMERS, NO_NEW_CLASHES;
extern int   SKIP_LOCAL_CONTACT_RANGE, SKIP_BB_CONTACT_RANGE;
extern Float SIDECHAIN_NOISE;
extern short READ_POTENTIAL, READ_DENSITY, USE_GO_POTENTIAL;

/* Alignment Variables */
struct align_cutoff {
    Float a;
    Float b;
};
extern struct align_cutoff **align_con_dist;
extern Float               **align_hard_core;
extern struct atom          *struct_native;
extern struct residue       *struct_residue;
extern int                   struct_natoms, struct_nresidues;
extern short                 struct_ncontacts, struct_nclashes;
extern int                  *seq_to_struct, *struct_to_seq;
extern int                  *map_to_seq, *map_to_struct;
extern struct segment {
    int a;
    int b;
} str_segment[100], seq_segment[100];
extern int  nseg;
extern int *helix;

/* Loop Variables */
extern int n_movable_residues;
extern int movable_residue_map[MAX_RES];  // This maps the movable residues to an array of size
                                          // n_movable_residues
extern int    res_atomno[MAX_RES];
extern short  yang_rotated_natoms;
extern short *yang_rotated_atoms;
extern char  *yang_not_rotated;
extern double a_PCA[MAXSEQUENCE], a_bCA[MAXSEQUENCE], phim[MAXSEQUENCE], psim[MAXSEQUENCE];
extern short  torsion_E[MAXSEQUENCE][6][6][6][6];
extern short  aromatic_E[9];
extern int    Naromatic;
extern int    Res_aromatic[MAXSEQUENCE];

/* Roatation variables */
extern Float rot_mat_00, rot_mat_10, rot_mat_20, rot_mat_01, rot_mat_11, rot_mat_21, rot_mat_02,
    rot_mat_12, rot_mat_22;

/* Debugging Variables */
#if DEBUG
extern unsigned char **debug_contacts;
extern unsigned char **debug_dcontacts;
extern short           debug_ncontacts;
extern unsigned char **debug_clashes;
extern short           debug_nclashes;
#endif

// pdb utils variable
extern int first_atom_res;

// fold variables
extern Float      native_rms, rms_RMSDmin, rms_Emin;
extern long int   mcstep_Emin, mcstep_RMSDmin;
extern float      delta_E, delta_T, delta_N, delta_all;
extern MPI_Status mpi_status;
extern float      E_RMSDmin;
extern float      E_RMSDmin_pot, E_RMSDmin_hbond, E_RMSDmin_tor, E_RMSDmin_sct, E_RMSDmin_aro,
    E_RMSDmin_constraint;

/*backbone variables*/
extern char path_dir[500];

/* Constraint Variables */
extern int    Max_N_constraints, N_constraints;
extern int  **constraint_array;
extern float *distances, *constraint_weights, k_constraint, E_constraint, prev_E_constraint,
    dE_constraint, mean_constraint_distance;

/* MPI Buffers */
extern float *Enode, *Tnode, *buf_in, *buf_out;
extern int   *Nnode, *Cnode;

#endif