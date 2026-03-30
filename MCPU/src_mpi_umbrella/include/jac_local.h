#ifndef JAC_LOCAL_H
#define JAC_LOCAL_H

// subroutine loop_Jacobian(r_n, r_ca, r_c, Jac)
void   loop_Jacobian(double r_n[3][3], double r_ca[3][3], double r_c[3][3], double *Jac);
double det3(double j1[3], double j2[3], double j3[3]);

#endif