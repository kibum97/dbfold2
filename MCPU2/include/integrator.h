#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include "topology.h"

class Integrator {
    public:
    Integrator();
    ~Integrator();

    Eigen::Matrix3Xd partialRotation(Eigen::Matrix3Xd& positions, const std::vector<size_t>& atomIDs, double angle, Eigen::Vector3d axis);
    void concertedRotation(Eigen::Matrix3Xd& positions, const std::vector<size_t>& atomIDs, double angle, Eigen::Vector3d axis);
    
    private:
    std::vector<size_t> selectAtomIDs(size_t atomID1, size_t atomID2, const Topology& topology);
    Eigen::Matrix3d createRotationMatrix(double angle, Eigen::Vector3d axis);
    Eigen::Quaterniond createQuaternion(double angle, Eigen::Vector3d axis);
};

#endif // INTEGRATOR_H