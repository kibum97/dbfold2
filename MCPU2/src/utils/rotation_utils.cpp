#include "utils/rotation_utils.h"

Eigen::AngleAxisf create_rotation_matrix(float angle, Eigen::Vector3f axis) {
    return Eigen::AngleAxisf(angle, axis);
}

Eigen::Quaternionf create_rotation_quaternion(float angle, Eigen::Vector3f axis) {
    return Eigen::Quaternionf(Eigen::AngleAxisf(angle, axis));
}

Eigen::Vector3f rotate_vector(const Eigen::Vector3f& vector, const Eigen::Vector3f& origin, Eigen::AngleAxisf rotation) {
    return rotation * (vector - origin) + origin;
}

Eigen::Matrix3Xf rotate_matrix(const Eigen::Matrix3Xf& matrix, const Eigen::Vector3f& origin, Eigen::AngleAxisf rotation) {
    Eigen::Matrix3f rotmat = rotation.toRotationMatrix();
    return (rotmat * (matrix.colwise() - origin)).colwise() + origin;
}