#include "utils/geometry_utils.h"
#include <algorithm>
#include <iostream>
#include <cmath>

float compute_angle(const Eigen::Vector3f& a, const Eigen::Vector3f& b) {
    // a and b should be normalized vectors
    float dot = a.dot(b);
    dot = std::max(-1.0f, std::min(1.0f, dot));
    float angle = std::acos(dot);
    return angle;
}

float compute_length(const Eigen::Vector3f& a, const Eigen::Vector3f& b) {
    Eigen::Vector3f r = a - b;
    float length = r.norm();
    return length;
}

float compute_bond_angle(const Eigen::Vector3f& a, const Eigen::Vector3f& b, const Eigen::Vector3f& c) {
    Eigen::Vector3f r1 = a - b;
    Eigen::Vector3f r2 = c - b;

    r1 = r1.normalized();
    r2 = r2.normalized();
    
    // Compute the bond angle between two bonds
    float angle = compute_angle(r1, r2);
    return angle;
}

float compute_dihedral_angle(const Eigen::Vector3f& a, const Eigen::Vector3f& b, const Eigen::Vector3f& c, const Eigen::Vector3f& d) {
    // Compute the dihedral angle between four atoms
    Eigen::Vector3f b1 = a - b;
    Eigen::Vector3f b2 = c - b;
    Eigen::Vector3f b3 = c - d;
    Eigen::Vector3f c1 = b2.cross(b3);    
    Eigen::Vector3f c2 = b1.cross(b2);
    float p1 = b1.dot(c1);
    p1 *= b2.norm();
    float p2 = c1.dot(c2);
    float angle = std::atan2(p1, p2);
    return angle;
}

float compute_plane_angle(const Eigen::Vector3f& n1, const Eigen::Vector3f& ca1, const Eigen::Vector3f& o1,
                          const Eigen::Vector3f& n2, const Eigen::Vector3f& ca2, const Eigen::Vector3f& o2) {
    // Compute the plane angle
    Eigen::Vector3f normal1 = compute_plane_vector(n1, ca1, o1);
    Eigen::Vector3f normal2 = compute_plane_vector(n2, ca2, o2);
    float angle = compute_angle(normal1, normal2);
    return angle;
}

float compute_bisect_angle(const Eigen::Vector3f& n1, const Eigen::Vector3f& ca1, const Eigen::Vector3f& o1,
                         const Eigen::Vector3f& n2, const Eigen::Vector3f& ca2, const Eigen::Vector3f& o2) {
    // Compute the bisect angle
    Eigen::Vector3f bisect1 = compute_bisect_vector(n1, ca1, o1);
    Eigen::Vector3f bisect2 = compute_bisect_vector(n2, ca2, o2);
    float angle = compute_angle(bisect1, bisect2);
    return angle;
}

float compute_delta_dihedral_angle(Eigen::Vector3f& a, Eigen::Vector3f& b, Eigen::Vector3f& c, Eigen::Vector3f& d_old, Eigen::Vector3f& d_new) {
    // Compute the delta dihedral angle between four atoms
    float angle_old = compute_dihedral_angle(a, b, c, d_old);
    float angle_new = compute_dihedral_angle(a, b, c, d_new);
    float delta_angle = angle_new - angle_old;
    return delta_angle;
}

Eigen::Vector3f compute_center(const std::vector<Eigen::Vector3f>& atoms) {
    Eigen::Vector3f center(0, 0, 0);
    for (const auto& atom : atoms) {
        center += atom;
    }
    center /= atoms.size();
    return center;
}

Eigen::Vector3f compute_plane_vector(const Eigen::Vector3f& a, const Eigen::Vector3f& b, const Eigen::Vector3f& c) {
    // Compute the plane vector between three atoms
    Eigen::Vector3f ba = b - a;
    Eigen::Vector3f bc = b - c;
    Eigen::Vector3f normal = ba.cross(bc);
    return normal.normalized();
}

Eigen::Vector3f compute_bisect_vector(const Eigen::Vector3f& a, const Eigen::Vector3f& b, const Eigen::Vector3f& c) {
    // Compute the bisect vector between three atoms
    Eigen::Vector3f ba = b - a;
    Eigen::Vector3f bc = b - c;
    ba = ba.normalized();
    bc = bc.normalized();
    Eigen::Vector3f bisect = (ba + bc) / 2.0f;
    return bisect.normalized();
}

float normalize_angle(float angle) {
    while (angle < 0) {
        angle += 360.0;
    }
    while (angle >= 360.0) {
        angle -= 360.0;
    }
    return angle;
}

int angle2index(float angle, float bin_size) {
    // Convert Radian to Degree
    angle = angle * (180.0f / M_PI);
    angle += 180.0f;
    return static_cast<int>(std::fmod(angle + 15.0f, 360.0f)/bin_size);
}