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

float prob_ang[20][100];

// int max_N_cysteines=100; //Ideally would figure out the number of cysteines by getting the number
// of unique residues in the cons traints list..for now assume no more than 100 /* 300>500 */

/* VZ:
 * The following two used to be local vars in fold.h, and struct_f2 also local to move.h
 * struct_f1 holds reference structure backbone.
 * I've moved them to here because I'd like the functions in move.h to access them
 * Yes even more global namespace pollution unfortunately.
 * A struct backbone holds five struct vector's named N, CA, C, CB, O
 * MAXSEQUENCE is macro for 350
 */

// struct backbone constraint_struct_f1[MAXSEQUENCE], constraint_struct_f2[MAXSEQUENCE];
Float Rg; /* radius of gyration */
Float sc_rms, new_rms;

/* globals for move.h */

short total_pairs, total_pairs2;
short TOTAL_SINGLE_LOOP_MOVES, TOTAL_DOUBLE_LOOP_MOVES, TOTAL_TRIPLE_LOOP_MOVES;
float jacobi_before, jacobi_after;
int   soln_no_before, soln_no_after;
int   total_ntorsions;

/* ----------------------- */

/*Bias*/
double k_bias;

short **substructures;  // AB

int umbrella;               /*Is umbrella simulation turned on?*/
int number_of_contacts_max; /*Number of contacts setpoint for first node*/
int contacts_step;          /*By how much does stepoint decrease for each node? */

int new_natives;

int  n_substructures;
int *substructure_sizes;  // AB

Float  dE_sct, prev_E_sct, E_sct;
Float  dE_aro, prev_E_aro, E_aro;
Float  dE_Rg, prev_E_Rg, E_Rg;
double cur_phi[MAXSEQUENCE], cur_psi[MAXSEQUENCE];

float cluster_phi[20][NOCLUSTERS];
float cluster_psi[20][NOCLUSTERS];

/* parameter file names */

/* vzhao: I increased the length limits to 500 from 100/150 for these: */

char substructure_path[500];  // AB
/* file names added 12DEC02 IAH */
char min_etot_file[500], min_drms_file[500];

/*Things defined by AB related to knowledge-based moves*/

/*Added by AB to allow myrank to be offset, for instnace when running multi-core unfolding
 * simulations without MPI...by default, it is not*/
int my_rank_offset = 0;

/* includes */
#include "align.h"
#include "constraint.h" /*AB on 11/26/19*/
#include "contacts.h"
#include "debug.h"
#include "energy.h"
#include "fold.h"
#include "globals.h"
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
