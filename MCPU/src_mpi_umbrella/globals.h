#ifndef GLOBALS_H
#define GLOBALS_H

#include "define.h"   /* for definitions */
#define pi 3.141592653589793238462643383279502884197e0
#define rad2deg 180.0e0 / pi

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>      // For strcmp

#include "misc_util.h" /* for Tiny */
#include "lattice_util.h" /* for struct cell */
#include "vector.h" /* for struct vector */

// Global Variables (extern tells the compiler they exist elsewhere)
extern int natoms;
extern FILE *STATUS;
extern int MAX_TYPES;
extern char PROTEIN_NAME[100]; // or whatever the size was
extern double mu;
extern double **potential;

/* These depend on your struct definitions, so make sure 
   the structs (atom, contact) are defined above these lines */
extern struct atom *native;
extern struct contact_data {
    unsigned char contacts : 1;
    unsigned char delta_contacts : 1;
    unsigned char clashes : 1;
    unsigned char delta_clashes : 1;
    unsigned char disulfide : 1;
    unsigned char check_contacts : 1;
    unsigned char check_clashes : 1;
    unsigned char check_hbond : 1;
    long int      hbond;
    long int      delta_hbond;
    int           closehb;  // if 0, r<2.5
    int           delta_closehb;
} **data, **struct_data;

struct int_vector {
    long int x, y, z;
};

struct vector {
    double x, y, z;
};


extern struct amino {
    char  name[4];
    char  symbol[2];
    char  torsion[4][4][4];
    Float avg_angle[4][4][4][4][4];
    int   ntorsions;
    int   nrotamers;
    int   rotate_natoms[4];
    char  rotate_atom[4][10][4];
};

struct atom {
    struct vector     xyz;
    struct int_vector xyz_int;
    char              res[5];  /* residue name */
    short             res_num; /* residue number */
    char              atomname[5];
    short             atomtype;
    short             smogtype;
    short             is_core;
    short             is_designed;
    short             is_sidechain;
    short             sec_structure;
    unsigned char     X, Y, Z;
    struct cell      *matrix;
};

struct residue {
    char    res[4];
    int     amino_num;
    short   is_core;
    short   is_designed;
    Float   psi;
    Float   phi;
    int     ntorsions;
    int     nrotamers;
    Float   avg_angle[4][4][4][4][4];
    short **rot_position;
    Float   chi[4];
    Float   native_chi[4];
    Float   tmpchi[4];
    int     CA, N, C, O, CB, CG, CE1, CE2, CZ2, CZ3;
    int     atomnumber[40];
};
struct residue *native_residue;

extern double rot_mat[3][3];

/* contact matrices */
extern short nclashes, ncontacts;
extern short delta_nclashes, delta_contacts;


/* globals for contacts.h */
extern long int          distance;
extern long int          X_int, Y_int, Z_int;
extern short             M, N, O, P;
extern short            *A, *B;
extern struct int_vector temp_xyz_int;
extern struct atom      *temp_atom, *temp_prev_atom;
extern struct cell      *temp_cell, *temp_cell3;
extern struct cell     **temp_cell_array;
extern unsigned char    *is_rotated;
extern short                **type_contacts;
extern long int **hard_core;
extern struct cutoff {
    long int a;
    long int b;
}            **contact_distance;

extern struct monte_carlo_flags {
    unsigned char clashed : 1;
    unsigned char init : 1;
} mc_flags;

extern int            MAX_TYPES, natom_type_list, bb_O_type, bb_N_type, bb_OXT_type;


struct pair {
    short a;
    short b;
};

struct pair *ab, *cd;

/* protein data */
/* These will be pointers to an array of struct atoms, representing the protein atoms */
extern struct atom    *native, *prev_native, *orig_native, *native_Emin, *native_RMSDmin;

extern int   natoms, nresidues; /* The number of atoms and residues in the structure */

extern int     min_seq_sep; /*How far apart in sequence do two residues have to be to qualify as a contact?
                        Only relevant for native contact calculation*/


/* MC variables */
extern long int mcstep, mcstep_Emin, mcstep_RMSDmin;
extern int      sidechain_step;
extern int     *cur_rotamers, old_rotamer;
extern Float    low_E, prev_E, E, native_E, Emin, Emin_pot, Emin_hbond, Emin_tor, Emin_sct, Emin_aro,
    Emin_constraint;
extern float E_RMSDmin_pot, E_RMSDmin_hbond, E_RMSDmin_tor, E_RMSDmin_sct, E_RMSDmin_aro,
    E_RMSDmin_constraint;
extern float E_RMSDmin;
extern Float dE, dE_pot, prev_E_pot, E_pot;
extern int   naccepted, n_sidechain_accepted;
extern int   nrejected, nothers, nomove;
extern int   sidemovedone;
struct monte_carlo {
    int   backbone_move;
    int   sel_res_num;
    int   sel_rotamer;
    int   is_phi;
    int   torsion;
    Float delta_angle[4];
    Float delta_phi_angle[3];
    Float delta_psi_angle[3];
    int   loop_size;
    int   sel_triplet;
    int   selected[3];
} mc;

/* Constraint variables--by AB!!! */
int   Max_N_constraints = 1000;
int   N_constraints;
int **constraint_array;  // We will declare a pointer to an array of Max_N_constraints pointers to
                         // ints...the ** indicates pointer to a pointer

float *distances;  // A pointer to an array of N_constraints distances

float *constraint_weights;       // A pointer to an array
int   *disulfide_pairs_attempt;  // The nth element will tell you whether the pair of residues
                               // indicated by the nth constraint is involved in a disulfide bond in
                               // the proposed configuration
extern int *disulfide_pairs;  // The nth element will tell you whether the pair of residues indicated by
                       // the nth constraint is involved in a disulfide bond in the currently
                       // accepted configuration
extern char  constraint_file[500];
extern char  rmsd_constraint_file[500];
extern char  movable_region_file[500];  // KP added for movable region constraint
float k_constraint              = 0;
float E_constraint;
float prev_E_constraint;
float dE_constraint;
// float E_constraint_now;
float mean_constraint_distance = 0;

extern int     natives;                /*Current number of native contacts for a given node*/

/* Sequencing and Alignment */
extern int nalign;
extern short *map_to_seq;
extern struct backbone *struct_f1, *struct_f2;

/* Energy and Math Counters */
extern int total_ntorsions, total_pairs, total_pairs2, total_hbond_pairs;
extern double new_rms, dE_hbond, dE_tor, dE_aro, dE_sct;
extern double prev_E_tor, prev_E_sct, prev_E_aro, prev_E_hbond;
extern double weight_hbond, weight_potential, weight_clash;
extern double MC_TEMP, k_bias;
extern int soln_no_before, jacobi_after, jacobi_before;

/* Rotation and Move Arrays */
extern int **rotate_natoms;
extern short ***rotate_atom;
extern char ***not_rotated;
extern int **loop_rotate_natoms;
extern short ***loop_rotate_atoms;
extern char ***loop_not_rotated;
extern int **rotate_sidechain_natoms;
extern short ***rotate_sidechain_atom;
extern char ***sidechain_not_rotated;

/* Constants and States */
extern int SIDECHAIN_MOVES;
extern int TOTAL_TRIPLE_LOOP_MOVES, TOTAL_SINGLE_LOOP_MOVES, TOTAL_DOUBLE_LOOP_MOVES;
extern double deg2rad;
extern char secstr[];

extern struct triplet {
    short a;
    short b;
    short c;
};

extern int new_natives;
extern short  *orig_contactstring;
extern int all_rotated_natoms;
extern short *all_rotated_atoms;
extern int move_cycle;
extern struct triplet *residue_triplets;
extern double **cluster_phi;
extern double **cluster_psi;
extern double *cur_phi;
extern double *cur_psi;

extern struct angles {
    Float chis[100][4];
};
extern int *no_chi_list;
extern double **prob_ang;
extern int ***sidechain_torsion;
extern struct angles *rotamer_angles;
extern double ***deviation_ang;

extern struct atom_type {
    char res_name[4];
    char atom_name[4];
    int  type_num;
};

#endif