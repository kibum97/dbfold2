#ifndef STURM_H
#define STURM_H

#include "globals.h"

/* Modified version of the code downloaded
 * from http://www.acm.org/pubs/tog/GraphicsGems/gems/Sturm/
 * Modified by Chaok Seok 2003.
 *
 * Using Sturm Sequences to Bracket Real Roots of Polynomial Equations
 * by D.G. Hook and P.R. McAree
 * from "Graphics Gems", Academic Press, 1990
 */

#ifdef __cplusplus
extern "C" {
#endif

#define PRINT_LEVEL 0
#define MAX_ORDER 16
#define MAXPOW 32
#define SMALL_ENOUGH 1.0e-18

/*
 * structure type for representing a polynomial
 */
typedef struct p {
    int    ord;
    double coef[MAX_ORDER + 1];
} poly;

extern double RELERROR;
extern int    MAXIT, MAX_ITER_SECANT;

void       initialize_sturm(double *tol_secant, int *max_iter_sturm, int *max_iter_secant);
void       solve_sturm(int *p_order, int *n_root, double *poly_coeffs, double *roots);
double     hyper_tan(double a, double x);
static int modp(poly *u, poly *v, poly *r);
int        buildsturm(int ord, poly *sseq);
int        numroots(int np, poly *sseq, int *atneg, int *atpos);
int        numchanges(int np, poly *sseq, double a);
void       sbisect(int np, poly *sseq, double min, double max, int atmin, int atmax, double *roots);
double     evalpoly(int ord, double *coef, double x);
int        modrf(int ord, double *coef, double a, double b, double *val);

#ifdef __cplusplus
}
#endif

#endif
