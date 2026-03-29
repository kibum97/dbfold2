#ifndef ROTATION_UTILS_H
#define ROTATION_UTILS_H

#include <Eigen/Dense>

Eigen::AngleAxisf  create_rotation_matrix(float angle, Eigen::Vector3f axis);
Eigen::Quaternionf create_rotation_quaternion(float angle, Eigen::Vector3f axis);
Eigen::Vector3f    rotate_vector(const Eigen::Vector3f &vector, const Eigen::Vector3f &origin,
                                 Eigen::AngleAxisf rotation);
Eigen::Matrix3Xf   rotate_matrix(const Eigen::Matrix3Xf &matrix, const Eigen::Vector3f &origin,
                                 Eigen::AngleAxisf rotation);

#endif  // ROTATION_UTILS_H