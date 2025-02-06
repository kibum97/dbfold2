#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include "topology.h"
#include "system.h"


class Integrator {
    public:
    Integrator();
    ~Integrator();

    Eigen::Matrix3Xd partialRotation(Eigen::Matrix3Xd& positions, const std::vector<size_t>& atomIDs, double angle, Eigen::Vector3d axis);
    Eigen::Matrix3Xd sidechainRotation(size_t resID, Eigen::Matrix3Xd& positions, std::unordered_map<size_t, std::vector<RotamerData>> rotamerIDMap, std::unordered_map<size_t, std::vector<TorsionIDData>> rotatingSCAtomsMap, std::vector<std::array<double, 4>> sidechain_torsions);
    void concertedRotation(Eigen::Matrix3Xd& positions, const std::vector<size_t>& atomIDs, double angle, Eigen::Vector3d axis);
    private:
    std::vector<size_t> selectAtomIDs(size_t atomID1, size_t atomID2, const Topology& topology);
    Eigen::Matrix3d createRotationMatrix(double angle, Eigen::Vector3d axis);
    Eigen::Quaterniond createQuaternion(double angle, Eigen::Vector3d axis);
};

#endif // INTEGRATOR_H