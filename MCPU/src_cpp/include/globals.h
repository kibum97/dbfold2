#ifndef GLOBALS_H
#define GLOBALS_H

#include <Eigen/Dense>

namespace MCPU2 {
    #ifdef USE_SINGLE_PRECISION
        using Real = float;
    #else
        using Real = double;
    #endif

    using Vec3  = Eigen::Vector3<Real>;
    using Vec3i = Eigen::Vector3i;
    using Mat3  = Eigen::Matrix3<Real>;

    // One-letter cast helper to keep math code from getting cluttered
    template <typename T>
    inline auto asReal(const T& eigen_obj) { return eigen_obj.template cast<Real>(); }
}

using namespace MCPU2;
#define _USE_MATH_DEFINES
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "define.h"

/* Math constants */
static const double pi      = M_PI;
static const double rad2deg = 180.0 / M_PI;
static const double deg2rad = M_PI / 180.0;

/* Struct definitions */
struct pair {
    short a, b;
};
struct cutoff {
    long int a, b;
};
struct triplet {
    short a, b, c;
};
struct angles {
    float chis[100][4];
};
struct align_cutoff {
    float a, b;
};
struct segment {
    int a, b;
};

struct atom {
    Vec3     xyz;
    Vec3i    xyz_int;
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
    float avg_angle[4][4][4][4][4];
    int   ntorsions;
    int   nrotamers;
    int   rotate_natoms[4];
    char  rotate_atom[4][10][4];
};

struct residue {
    char    res[4];
    int     amino_num, ntorsions, nrotamers;
    short   is_core, is_designed;
    float   psi, phi, chi[4], native_chi[4], tmpchi[4];
    float   avg_angle[4][4][4][4][4];
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
    Vec3 N;
    Vec3 CA;
    Vec3 C;
    Vec3 CB;
    Vec3 O;
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
    float delta_angle[4], delta_phi_angle[3], delta_psi_angle[3];
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

struct fragments {
    int x1;
    int x2;
};

struct alignment {
    struct fragments seqptr[MAXFRAG];
    struct fragments structptr[MAXFRAG];
    int              NFRAG;
};


/* Topology - Static biological structure */
struct Topology {
    int natoms, nresidues;
    int first_atom_res;
    int natom_type_list;
    int total_ntorsions;
    int Naromatic;
    int n_movable_residues;
    int bb_O_type, bb_N_type, bb_OXT_type;
    int min_seq_sep;
    int nalign;

    /* Alignment Targets */
    int struct_natoms, struct_nresidues;

    char secstr[MAXSEQUENCE];
    int  *map_to_seq, *map_to_struct;
    int  *seq_to_struct, *struct_to_seq;
    int  no_chi_list[20];
    int  *helix;
    int  Res_aromatic[MAXSEQUENCE];
    int  movable_residue_map[MAX_RES];
    int  is_template[MAX_RES];
    int  res_atomno[MAX_RES];

    struct atom_type atom_type_list[200];
    struct amino    *amino_acids;
};

/* System - Physics weights and potential parameters */
struct System {
    int MAX_TYPES;
    int USE_GO_POTENTIAL;
    int SEQ_DEP_HB;

    float    weight_potential, weight_clash, weight_rms, weight_hbond;
    float    NATIVE_ATTRACTION, NON_NATIVE_REPULSION, NON_SPECIFIC_ENERGY;
    float    ALPHA, LAMBDA;
    float    mu, beta;
    double   k_bias;
    float    rmsd_constraint;
    double   contact_calpha_cutoff;
    long int NO_HBOND;

    /* Tables and Matrices (Future Eigen Targets) */
    float          *radii;
    float         **potential;
    long int      **hard_core;
    struct cutoff **contact_distance;

    double **prob_ang;
    double **cluster_phi;
    double **cluster_psi;
    float    deviation_ang[20][100][4];

    short      torsion_E[MAXSEQUENCE][6][6][6][6];
    short      aromatic_E[9];
    short *****sct_E;
    float      seq_hb[3][20][20];

    /* Constraints Definitions */
    int    Max_N_constraints, N_constraints;
    int  **constraint_array;
    float *distances;
    float *constraint_weights;
    float  k_constraint;

    /* Alignment Restraints */
    struct align_cutoff **align_con_dist;
    float               **align_hard_core;
    struct segment        str_segment[100], seq_segment[100];
    int                   nseg;
};

/* Monte Carlo Integrator - How the system evolves */
struct MCIntegrator {
    float MC_TEMP;
    float MC_TEMP_MIN;
    float TEMP_STEP;
    float STEP_SIZE;
    float YANG_SCALE;
    float CLUSTER_MOVE;
    float SIDECHAIN_NOISE;

    float YANG_MOVE;
    float USE_CLUSTER;
    int   USE_ROTAMERS;
    int   USE_ROT_PROB;
    int   USE_SIDECHAINS;
    float USE_GLOBAL_BB_MOVES;
    int   SIDECHAIN_MOVES;
    int   MAX_CLUSTERSTEP;

    float LATTICE_SIZE;
    int   MATRIX_SIZE;
    int   HALF_MATRIX_SIZE;

    int TOTAL_TRIPLE_LOOP_MOVES;
    int TOTAL_SINGLE_LOOP_MOVES;
    int TOTAL_DOUBLE_LOOP_MOVES;
    int total_triplets;

    struct triplet *residue_triplets;
    char         ***not_rotated;
    short        ***rotate_atom;
    short         **rotate_natoms;
    short        ***loop_rotate_atoms;
    short         **loop_rotate_natoms;
    char         ***loop_not_rotated;
    short        ***sidechain_torsion;
    struct angles  *rotamer_angles;
    short         **rotate_sidechain_natoms;
    short        ***rotate_sidechain_atom;
    char         ***sidechain_not_rotated;
    short         **loop_int_rotate_natoms;
    short        ***loop_int_rotate_atoms;

    short  yang_rotated_natoms;
    short *yang_rotated_atoms;
    char  *yang_not_rotated;
};

/* Context - The living state of the simualtion */
struct Context {
    long int mcstep;
    int      move_cycle;
    int      sel_num;
    int      soln_no_before, jacobi_after, jacobi_before;

    /* Coordinates and State Flags */
    struct atom             *native, *prev_native, *orig_native, *native_Emin, *native_RMSDmin;
    struct residue          *native_residue;
    struct atom             *struct_native;
    struct residue          *struct_residue;
    struct backbone          struct_f1[MAXSEQUENCE], struct_f2[MAXSEQUENCE];
    struct monte_carlo       mc;
    struct monte_carlo_flags mc_flags;

    /* Dynamic Arrays */
    short **type_contacts;
    int    *cur_rotamers;
    short  *all_rotated_atoms;
    int    *disulfide_pairs;
    int    *disulfide_pairs_attempt;
    short  *orig_contactstring;
    double *cur_phi;
    double *cur_psi;
    double  a_PCA[MAXSEQUENCE], a_bCA[MAXSEQUENCE], phim[MAXSEQUENCE], psim[MAXSEQUENCE];

    /* Trackers */
    short nclashes, ncontacts, delta_nclashes, delta_contacts;
    short struct_ncontacts, struct_nclashes;
    int   total_pairs, total_pairs2, total_hbond_pairs;
    int   diff_number_of_contacts_new;
    int   diff_number_of_contacts_current;
    int   new_natives;
    int   natives;
    int   sidechain_step;
    int   backbone_accepted;
    int   old_rotamer;
    int   sidemovedone;
    short all_rotated_natoms;

    /* Current Matrix and Pairs */
    struct contact_data  **data, **struct_data;
    struct cell         ***the_matrix;
    struct pair           *ab, *cd;
    struct pair           *hbond_pair;
    struct hydrogen_bond **hbonds;

    /* Current Energies */
    float  native_E;
    float  E, E_pot, dE, dE_pot, prev_E, prev_E_pot;
    float  E_tor, E_sct, E_aro, E_Rg, dE_tor, dE_sct, dE_aro, dE_Rg;
    float  prev_E_aro, prev_E_sct, prev_E_tor;
    float  prev_E_hbond, E_hbond, dE_hbond;
    float  hydrogen_bond;
    short *hbond_E;
    int    helix_sheet;
    float  E_constraint, prev_E_constraint, dE_constraint, mean_constraint_distance;
    float  delta_E, delta_T, delta_N, delta_all;

    /* Minimums tracking */
    float Emin, Emin_pot, Emin_hbond, Emin_tor, Emin_sct, Emin_aro, Emin_constraint;
    float E_RMSDmin, E_RMSDmin_pot, E_RMSDmin_hbond, E_RMSDmin_tor, E_RMSDmin_sct, E_RMSDmin_aro,
        E_RMSDmin_constraint;
    float    rms_RMSDmin, rms_Emin, native_rms;
    long int mcstep_Emin, mcstep_RMSDmin;
};

/* Simulation - The wrapper, I/O and Replica Exchange */
struct Simulation {
    long int seed;
    long int MC_STEPS;
    long int MC_INIT_STEPS;
    long int MC_ANNEAL_STEPS;
    long int MC_PRINT_STEPS;
    long int MC_PDB_PRINT_STEPS;
    long int MC_REPLICA_STEPS;

    /* Options */
    int   PRINT_PDB;
    short READ_POTENTIAL, READ_DENSITY;
    int   NO_NEW_CLASHES;
    int   SKIP_LOCAL_CONTACT_RANGE, SKIP_BB_CONTACT_RANGE;
    int   umbrella, number_of_contacts_max;
    int   number_of_contacts_setpoint;
    int   contacts_step;

    /* MPI / Replica Exchange */
    int        nprocs, myrank, ierr;
    MPI_Comm   mpi_world_comm;
    MPI_Status mpi_status;
    int        NODES_PER_TEMP;
    int        my_rank_offset;
    int        current_replica;
    int        MAX_EXCHANGE;
    int       *accepted_replica, *rejected_replica, *attempted_replica;
    int       *replica_index;
    float     *Enode, *Tnode;
    Eigen::Matrix3Xd buf_in, buf_out;
    int       *Nnode, *Cnode;

    /* Statistics */
    int naccepted, n_sidechain_accepted;
    int nrejected, nothers, nomove;

    /* File handles and paths */
    FILE *fpdb;
    FILE *STATUS;
    FILE *DATA;
    char  path_dir[1024];
    char  native_file[1024], structure_file[1024], triplet_file[1024];
    char  sctorsion_file[1024], sec_str_file[1024], template_file[1024];
    char amino_data_file[1024], rotamer_data_file[1024], potential_file[1024], atom_type_file[1024];
    char helicity_data[1024], hydrogen_bonding_data[1024];
    char substructure_path[1024], alignment_file[1024];
    char std_prefix[1024];
    char aromatic_file[1024];
    char PROTEIN_NAME[100];
    char pdb_out_file[1024];
    char std_file[1024];
    char native_directory[1024];
    char constraint_file[1024];
    char rmsd_constraint_file[1024];
    char movable_region_file[1024];

    /* Debug */
#if DEBUG
    unsigned char **debug_contacts;
    unsigned char **debug_dcontacts;
    unsigned char **debug_clashes;
    short           debug_ncontacts;
    short           debug_nclashes;
#endif
};

extern struct Topology     top;
extern struct System       sys;
extern struct MCIntegrator integrator;
extern struct Context      ctx;
extern struct Simulation   sim;

extern float             rot_mat[3][3];
extern struct atom      *temp_atom, *temp_prev_atom;
extern Vec3i             temp_xyz_int;
extern struct cell      *temp_cell, *temp_cell3;
extern struct cell     **temp_cell_array;
extern unsigned char    *is_rotated;
extern long int          distance;

extern float  d_memory;
extern short  M, N, O, P;
extern short *A, *B;
extern float  rot_mat_00, rot_mat_10, rot_mat_20, rot_mat_01, rot_mat_11, rot_mat_21, rot_mat_02,
    rot_mat_12, rot_mat_22;

extern struct alignment constraint_align;

#endif

// extern float    hydrogen_bond;