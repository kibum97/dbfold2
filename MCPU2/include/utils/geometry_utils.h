#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <Eigen/Dense>

float compute_length(const Eigen::Vector3f &a, const Eigen::Vector3f &b);

float compute_angle(const Eigen::Vector3f &a, const Eigen::Vector3f &b);
float compute_bond_angle(const Eigen::Vector3f &a, const Eigen::Vector3f &b,
                         const Eigen::Vector3f &c);
float compute_dihedral_angle(const Eigen::Vector3f &a, const Eigen::Vector3f &b,
                             const Eigen::Vector3f &c, const Eigen::Vector3f &d);
float compute_plane_angle(const Eigen::Vector3f &n1, const Eigen::Vector3f &ca1,
                          const Eigen::Vector3f &o1, const Eigen::Vector3f &n2,
                          const Eigen::Vector3f &ca2, const Eigen::Vector3f &o2);
float compute_bisect_angle(const Eigen::Vector3f &n1, const Eigen::Vector3f &ca1,
                           const Eigen::Vector3f &o1, const Eigen::Vector3f &n2,
                           const Eigen::Vector3f &ca2, const Eigen::Vector3f &o2);
float compute_delta_dihedral_angle(Eigen::Vector3f &a, Eigen::Vector3f &b, Eigen::Vector3f &c,
                                   Eigen::Vector3f &d_old, Eigen::Vector3f &d_new);

Eigen::Vector3f compute_center(const std::vector<Eigen::Vector3f> &atoms);
Eigen::Vector3f compute_plane_vector(const Eigen::Vector3f &a, const Eigen::Vector3f &b,
                                     const Eigen::Vector3f &c);
Eigen::Vector3f compute_bisect_vector(const Eigen::Vector3f &a, const Eigen::Vector3f &b,
                                      const Eigen::Vector3f &c);
float           normalize_angle(float angle);
int             angle2index(float angle, float bin_size);

#endif  // GEOMETRY_H