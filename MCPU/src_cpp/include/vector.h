#ifndef VECTOR_H
#define VECTOR_H

#include <math.h>
#include "define.h"
#include "globals.h"

double D2(const Vec3& x, const Vec3& y); /* calculates the square of the distance */
void  Zero(Vec3* x);
void  MakeVector(const Vec3& x, const Vec3& y, Vec3* z);
double Angle(Vec3& x, Vec3& y);
double Dot(const Vec3& x, const Vec3& y);
double Norm(const Vec3& x);
void  Copy(const Vec3& x, Vec3* y);
void  Inverse(Vec3* x);
void  Scale(double factor, Vec3* x);
void  Add(const Vec3& x, Vec3* y);
void  Normalize(Vec3* x);
void  CrossProduct(const Vec3& x, const Vec3& y, Vec3* z);
void  bisect(const Vec3& x, const Vec3& y, Vec3* z);
void  RotateX(double angle, Vec3* A);
void  RotateY(double angle, Vec3* A);
void  RotateZ(double angle, Vec3* A);
void  MakeRotationMatrix(double angle, const Vec3& C);
void  Rotate(const Mat3& rot_mat, Vec3* A);
double Angle(const Vec3& x, const Vec3& y);

#endif