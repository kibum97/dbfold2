#ifndef JAC_LOCAL_H
#define JAC_LOCAL_H

#include <array>
#include "globals.h"

inline constexpr int MAX_ATTEMPTS = 2;

void loop_Jacobian(Mat3 r_n, 
                   Mat3 r_ca, 
                   Mat3 r_c, 
                   double *Jac);

#endif