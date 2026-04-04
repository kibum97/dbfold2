#ifndef DEFINE_H
#define DEFINE_H

// NOTE THAT THERE ARE PARAMTERS IN GETRMS.H
// MAXFRAG & MAXSEQUENCE
#define MAX_ATOMS 8000
#define MAX_CELL_ATOMS 100
#define MAXSEQUENCE 1000
#define MAXFRAG 30

#define MAX_RES 500
#define MSAR 15  // Max Sidechain Atoms per Residue

typedef float Float;
#define PRECISION 1000000
#define INT_PRECISION 100
#define PI 3.141592

#ifndef DEBUG
#define DEBUG 0
#endif
#define COMPILE_FOLDING 1

#define NSMOGTYPES 13
#define DISTANCE_DEPENDENCE 0
#define POTENTIAL_BIN_SIZE 0.5
#define MAX_DIST_BINS 15
#define DIST_BIN_MULT 10000 /* 1 * INT_PRECISION * INT_PRECISION */
#define NUM_ATOMNAMES 36
#define MIN_LEN_TEMPLATE 6

#define POTNTL_WEIGHT 0.4
#define HBOND_WEIGHT 1.35
#define TOR_WEIGHT 1.35
#define SCT_WEIGHT 2.50
#define ARO_WEIGHT 5.0
#define AROMATIC_DISTANCE 7.0

// AB moved all the following to backbone.h and defined as proper variables...
// #define CLUSTER_MOVE		0.33  //originally 0.3333
// #define USE_CLUSTER		0.0       // 0.0 to turn knowledge based moves off, was 0.1
// originally (knowledge moves on) #define MAX_CLUSTERSTEP 0  //Added by AB...all MC steps after
// this will always have USE_CLUSTER set to 0
#define NOCLUSTERS 30
#define CLUSTER_NOISE 10.

#define CUT_SECSTR 4
#define RDTHREE_CON 2.0

#define beta_favor 3.0  // MODIFY THIS for FACTOR
#define HB_CUTOFF 2.5
#define HB_INNER 2.5
#define HB_PENALTY 1.0
#define ang_int 20.
#define ang3_int 20.
#define LANG 9
#define LTOR 9

#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "rng.h"

extern short three[5];
const static char *at_types[NSMOGTYPES] __attribute__((unused)) = {
    "CA", "C3", "C2", "CP", "CC", "OC", "OB", "OD", "NC", "ND", "NM", "SP", "ME"};

const static char *ATOMNAME[NUM_ATOMNAMES] __attribute__((unused)) = {
    "C",   "CA",  "CB",  "CD",  "CD1", "CD2", "CE",  "CE1", "CE2", "CE3", "CG",  "CG1",
    "CG2", "CH2", "CZ",  "CZ2", "CZ3", "N",   "ND1", "ND2", "NE",  "NE1", "NE2", "NZ",
    "O",   "OCT", "OD1", "OD2", "OE1", "OE2", "OG",  "OG1", "OH",  "OXT", "SD",  "SG"};

#endif