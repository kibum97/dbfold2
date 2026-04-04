#include "vector.h"
#include "misc_util.h"

double D2(const Vec3& x, const Vec3& y) {
    return (x - y).squaredNorm();
}

void Zero(Vec3* x) {
    x->setZero();
}

void MakeVector(const Vec3& x, const Vec3& y, Vec3* z) {
    *z = y - x;
}

double Angle(Vec3& x, Vec3& y) {
    // Re-using your Tiny() macro/function for numerical stability
    x.normalize();
    y.normalize();
    return std::acos(Tiny(x.dot(y)));
}

double Dot(const Vec3& x, const Vec3& y) { 
    return x.dot(y); 
}

double Norm(const Vec3& x) { 
    return x.norm(); 
}

void Copy(const Vec3& x, Vec3* y) {
    *y = x;
}

void Scale(double factor, Vec3* x) {
    *x *= factor;
}

void Inverse(Vec3* x) {
    *x = -(*x);
}

void Add(const Vec3& x, Vec3* y) {
    *y += x;
}

void Normalize(Vec3* x) {
    x->normalize();
}

void CrossProduct(const Vec3& x, const Vec3& y, Vec3* z) {
    *z = x.cross(y);
}

// C = (A + B) / 2
void bisect(const Vec3& x, const Vec3& y, Vec3* z) {
    *z = (x + y) / 2.0;
}

// 2D Rotation around Z axis
void RotateZ(double angle, Vec3* A) {
    *A = Eigen::AngleAxis<double>(angle, Vec3::UnitZ()) * (*A);
}

// 2D Rotation around X axis
void RotateX(double angle, Vec3* A) {
    *A = Eigen::AngleAxis<double>(angle, Vec3::UnitX()) * (*A);
}

// 2D Rotation around Y axis
void RotateY(double angle, Vec3* A) {
    *A = Eigen::AngleAxis<double>(angle, Vec3::UnitY()) * (*A);
}

// Creates a rotation matrix around axis C
void MakeRotationMatrix(double angle, const Vec3& C) {
    // Eigen's AngleAxis automatically handles the quaternion/axis-angle math 
    Mat3 rot_mat = Eigen::AngleAxis<double>(angle, C.normalized()).toRotationMatrix();
}

// Applies a rotation matrix to A
// If keeping the C-style array for rot_mat, use Eigen::Map.
// Assuming rot_mat is now a Mat3:
void Rotate(const Mat3& rot_mat, Vec3* A) {
    *A = rot_mat * (*A);
}
