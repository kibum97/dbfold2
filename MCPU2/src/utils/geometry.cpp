#include "utils/geometry.h"

double computeAngle(const Eigen::Vector3d& a, const Eigen::Vector3d& b) {
    Eigen::Vector3d normA = a.normalized();
    Eigen::Vector3d normB = b.normalized();
    double dot = normA.dot(normB);
    dot = std::max(-1.0, std::min(1.0, dot));
    double angle = std::acos(dot);
    return angle;
}

double computeDihedralAngle(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c, const Eigen::Vector3d& d) {
    // Compute the dihedral angle between four atoms
    Eigen::Vector3d b1 = a - b;
    Eigen::Vector3d b2 = c - b;
    Eigen::Vector3d b3 = c - d;
    Eigen::Vector3d c1 = b2.cross(b3);    
    Eigen::Vector3d c2 = b1.cross(b2);
    double p1 = b1.dot(c1);
    p1 *= b2.norm();
    double p2 = c1.dot(c2);
    double angle = std::atan2(p1, p2);
    return angle;
}

double computePlaneAngle(const Eigen::Vector3d& n1, const Eigen::Vector3d& ca1, const Eigen::Vector3d& o1,
                         const Eigen::Vector3d& n2, const Eigen::Vector3d& ca2, const Eigen::Vector3d& o2) {
    // Compute the plane angle between four atoms
    Eigen::Vector3d b1 = ca1 - n1;
    Eigen::Vector3d b2 = ca1 - o1;
    Eigen::Vector3d b3 = ca2 - n2;
    Eigen::Vector3d b4 = ca2 - o2;
    Eigen::Vector3d normal1 = b1.cross(b2);
    Eigen::Vector3d normal2 = b3.cross(b4);
    double angle = computeAngle(normal1, normal2);
    return angle;
}

double computeBisectAngle(const Eigen::Vector3d& n1, const Eigen::Vector3d& ca1, const Eigen::Vector3d& o1,
                         const Eigen::Vector3d& n2, const Eigen::Vector3d& ca2, const Eigen::Vector3d& o2) {
    // Compute the bisect angle between four atoms
    Eigen::Vector3d b1 = ca1 - n1;
    Eigen::Vector3d b2 = ca1 - o1;
    Eigen::Vector3d b3 = ca2 - n2;
    Eigen::Vector3d b4 = ca2 - o2;
    Eigen::Vector3d bisect1 = (b1 + b2) / 2.0;
    Eigen::Vector3d bisect2 = (b3 + b4) / 2.0;
    double angle = computeAngle(bisect1, bisect2);
    return angle;
}

Eigen::Vector3d computeCenter(const std::vector<Eigen::Vector3d>& atoms) {
    Eigen::Vector3d center(0, 0, 0);
    for (const auto& atom : atoms) {
        center += atom;
    }
    center /= atoms.size();
    return center;
}

double normalizeAngle(double angle) {
    while (angle < 0) {
        angle += 360.0;
    }
    while (angle >= 360.0) {
        angle -= 360.0;
    }
    return angle;
}