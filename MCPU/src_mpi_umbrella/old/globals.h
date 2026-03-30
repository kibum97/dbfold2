#ifndef GLOBALS_H
#define GLOBALS_H

#include "define.h"   /* for definitions */

/* Math constants */
#define pi 3.141592653589793238462643383279502884197e0
#define rad2deg 180.0e0 / pi

/* Include standard libraries */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>      // For strcmp
#include <mpi.h>       // For MPI functions

/* Include custom headers - may be replaced*/
#include "misc_util.h" /* for Tiny */
#include "lattice_util.h" /* for struct cell */

// Struct definitions and global variable declarations will go here.
// Make sure to include any necessary headers before this point,
// and to define any structs that are used in the global variable declarations.

// Vector related structs
struct int_vector {
    long int x, y, z;
};
struct vector {
    float x, y, z;
};
struct pair {
    short a;
    short b;
};

// Protein structure related structs
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



// Simulation/System related structs
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
};

struct monte_carlo_flags {
    unsigned char clashed : 1;
    unsigned char init : 1;
};

// Contact data struct
struct contact_data {
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
extern struct contact_data **data, **struct_data;


extern struct residue *native_residue;
extern struct monte_carlo mc;

extern Float rot_mat[3][3];

/* contact matrices */
extern short nclashes, ncontacts;
extern short delta_nclashes, delta_contacts;


/* globals for contacts.h */

extern long int          X_int, Y_int, Z_int;
extern short             M, N, O, P;
extern short            *A, *B;












extern Float         LATTICE_SIZE;
extern unsigned char MATRIX_SIZE, HALF_MATRIX_SIZE;
struct cell ***the_matrix;

struct pair *ab, *cd;

/* protein data */
/* These will be pointers to an array of struct atoms, representing the protein atoms */
extern struct atom    *native, *prev_native, *orig_native, *native_Emin, *native_RMSDmin;

extern int   natoms, nresidues; /* The number of atoms and residues in the structure */

 /*How far apart in sequence do two residues have to be to qualify as a contact?
                        Only relevant for native contact calculation*/


/* MC variables */
extern long int mcstep, 


extern Float    low_E, prev_E, E, native_E, Emin, Emin_pot, Emin_hbond, Emin_tor, Emin_sct, Emin_aro,
    Emin_constraint;


extern Float dE, dE_pot, prev_E_pot, E_pot;



/* Constraint variables--by AB!!! */
extern int   Max_N_constraints;
extern int   N_constraints;
extern int **constraint_array;  // We will declare a pointer to an array of Max_N_constraints pointers to
                         // ints...the ** indicates pointer to a pointer

extern float *distances;  // A pointer to an array of N_constraints distances

extern float *constraint_weights;       // A pointer to an array
  // The nth element will tell you whether the pair of residues
                               // indicated by the nth constraint is involved in a disulfide bond in
                               // the proposed configuration
  // The nth element will tell you whether the pair of residues indicated by
                       // the nth constraint is involved in a disulfide bond in the currently
                       // accepted configuration
  // KP added for movable region constraint
extern float k_constraint;
extern float E_constraint;
extern float prev_E_constraint;
extern float dE_constraint;
// float E_constraint_now;
extern float mean_constraint_distance;

                /*Current number of native contacts for a given node*/

/* Sequencing and Alignment */
extern int nalign;
Float  dE_tor, prev_E_tor, E_tor;
Float  dE_sct, prev_E_sct, E_sct;
Float  dE_aro, prev_E_aro, E_aro;
Float  dE_Rg, prev_E_Rg, E_Rg;



/* rotation data structures */
          /* number of atoms affected by the rotation */
            /* atoms affected by the rotation */
     /* number of atoms affected by the rotation */
      /* atoms affected by the rotation */
extern short  **loop_int_rotate_natoms; /* number of atoms affected by the rotation */
extern short ***loop_int_rotate_atoms;  /* atoms affected by the rotation */





 /* pointer to atoms rotated by the current step */

extern MPI_Comm   mpi_world_comm;

extern int  nprocs, myrank, sel_num, ierr;
extern int  current_replica; /* Added 2021-07-12: to track a set of coordinates across replica space */
extern int *accepted_replica, *rejected_replica, *replica_index;
extern int *attempted_replica; /* VZ added 2-18 */
// int **accepted_replica, **rejected_replica, *replica_index; //AB changed to 2D array
extern float *Enode, *Tnode;
extern int   *Nnode, *Cnode;

extern float *buf_in, *buf_out;

/* Constants and States */
extern int SIDECHAIN_MOVES;

extern char secstr[];



extern short *all_rotated_atoms;










extern int myrank;
extern int nprocs;
extern MPI_Comm mpi_world_comm;

/* static data structures */
extern Float     *radii;

extern struct pair           *hbond_pair;

// Move related variables
extern short three[5];





extern int MAX_CLUSTERSTEP;

extern int move_cycle;
extern short           total_triplets;
extern short             *****sct_E;

/* contact matrices for debugging */




/* hydrogen bonding */

extern short                **type_contacts;


extern int                    n_hbond_int;

extern int                    total_hbond_pairs;



extern int   SEQ_DEP_HB;






/* alignment data structures */








// loop variables






//main 
extern long int       seed;

#endif