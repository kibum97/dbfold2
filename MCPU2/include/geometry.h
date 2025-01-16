#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <Eigen/Dense>
#include <iostream>
#include <cmath>

double computeDihedralAngle(const Eigen::Vector3d& p1, const Eigen::Vector3d& p2, const Eigen::Vector3d& p3, const Eigen::Vector3d& p4) {
    // Define the vectors between the atoms
    Eigen::Vector3d b1 = p2 - p1;
    Eigen::Vector3d b2 = p3 - p2;
    Eigen::Vector3d b3 = p4 - p3;

    // Compute the normal vectors to the planes
    Eigen::Vector3d n1 = b1.cross(b2);
    Eigen::Vector3d n2 = b2.cross(b3);

    // Normalize the normal vectors
    n1.normalize();
    n2.normalize();

    // Compute the angle between the normal vectors
    double cosTheta = n1.dot(n2);
    double sinTheta = b2.normalized().dot(n1.cross(n2));

    // Compute the dihedral angle
    double theta = std::atan2(sinTheta, cosTheta);

    // Convert the angle from radians to degrees
    double thetaDegrees = theta * (180.0 / M_PI);

    return thetaDegrees;
}

Eigen::Vector3d bisectVectors(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2) {
    Eigen::Vector3d bisector = v1.normalized() + v2.normalized();
    bisector.normalize();
    return bisector;
}

Eigen::Vector3d bisectVectors(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2) {
    Eigen::Vector3d bisector = v1.normalized() + v2.normalized();
    bisector.normalize();
    return bisector;
}

#endif // GEOMETRY_H