#ifndef MISC_UTIL_H
#define MISC_UTIL_H

#include "define.h"

Float Tiny(Float x);
void  squeeze(char *, int);
Float GaussianNum();

Float Tiny(Float x) { return ((Float)((int)(x * PRECISION))) / PRECISION; }

#endif