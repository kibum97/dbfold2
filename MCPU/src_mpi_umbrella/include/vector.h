#ifndef VECTOR_H
#define VECTOR_H

#include <math.h>

#include "define.h"
#include "globals.h"

Float D2(struct vector, struct vector); /* calculates the square of the distance */
void  Zero(struct vector *);
void  MakeVector(struct vector, struct vector, struct vector *);
Float Dot(struct vector, struct vector);
Float Norm(struct vector);
void  Copy(struct vector, struct vector *);
void  Inverse(struct vector *);
void  Scale(Float, struct vector *);
void  Add(struct vector, struct vector *);
void  Normalize(struct vector *);
void  CrossProduct(struct vector, struct vector, struct vector *);
void  bisect(struct vector, struct vector, struct vector *);
void  RotateX(Float, struct vector *);
void  RotateY(Float, struct vector *);
void  RotateZ(Float, struct vector *);
void  MakeRotationMatrix(Float, struct vector);
void  Rotate(Float[3][3], struct vector *);
Float Angle(struct vector, struct vector);

#endif