#ifndef GETRMS_H
#define GETRMS_H

#include "define.h"
#include "globals.h"

float getrms(struct backbone struct1[MAXSEQUENCE], struct backbone struct2[MAXSEQUENCE], struct alignment algn);
float getrms_fast(struct backbone* struct1, struct backbone* struct2, 
                  struct alignment& algn, Eigen::Matrix3Xd& P, Eigen::Matrix3Xd& Q);

#endif