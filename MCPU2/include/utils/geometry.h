#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <Eigen/Dense>
#include <iostream>
#include <cmath>

double computeAngle(const Eigen::Vector3d& a, const Eigen::Vector3d& b);
double computeDihedralAngle(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c, const Eigen::Vector3d& d);
double computePlaneAngle(const Eigen::Vector3d& n1, const Eigen::Vector3d& ca1, const Eigen::Vector3d& o1,
                         const Eigen::Vector3d& n2, const Eigen::Vector3d& ca2, const Eigen::Vector3d& o2);
double computeBisectAngle(const Eigen::Vector3d& n1, const Eigen::Vector3d& ca1, const Eigen::Vector3d& o1,
                          const Eigen::Vector3d& n2, const Eigen::Vector3d& ca2, const Eigen::Vector3d& o2);
Eigen::Vector3d computeCenter(const std::vector<Eigen::Vector3d>& atoms);
double normalizeAngle(double angle);

#endif // GEOMETRY_H