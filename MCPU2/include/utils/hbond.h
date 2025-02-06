#ifndef HBOND_H
#define HBOND_H

#include <Eigen/Dense>
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include "utils/geometry.h"

int computeSecondaryStructure(size_t donorResID, std::array<size_t, 7> donorAtomIDs,
                              size_t acceptorResID, std::array<size_t, 8> acceptorAtomIDs,
                              Eigen::Matrix3Xd positions);
bool checkHydrogenBond(size_t donorResID, std::array<size_t, 7> donorAtomIDs,
                       size_t acceptorResID, std::array<size_t, 8> acceptorAtomIDs,
                       Eigen::Matrix3Xd positions);
std::array<double, 3> computeDistanceFeatures(size_t donorResID, std::array<size_t, 7> donorAtomIDs,
                                              size_t acceptorResID, std::array<size_t, 8> acceptorAtomIDs,
                                              Eigen::Matrix3Xd positions);
std::array<double, 4> computeDihedralFeatures(size_t donorResID, std::array<size_t, 7> donorAtomIDs,
                                              size_t acceptorResID, std::array<size_t, 8> acceptorAtomIDs,
                                              Eigen::Matrix3Xd positions);
Eigen::Vector3d inferHydrogenPosition(const Eigen::Vector3d& c, const Eigen::Vector3d& n, const Eigen::Vector3d& ca);
int computeHBondindexfromAngles(int ss, std::array<double, 3> bisect_angles, std::array<double, 3> plane_angles);
int computeHBondindex(int ss, std::array<int, 3> bisect_indices, std::array<int, 3> plane_indices);

#endif // HBOND_H