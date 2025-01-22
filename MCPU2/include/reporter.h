#ifndef REPORTER_H
#define REPORTER_H

#include <iostream>
#include <Eigen/Dense>
#include "xtc/xdrfile.h"
#include "xtc/xdrfile_xtc.h"

class Reporter {
public:
    Reporter(const char* filename, int natoms, float timestep, float precision, const Eigen::Matrix3d& periodicBox);
    ~Reporter();
    bool writeTrajectory(int frame_idx, Eigen::Matrix3Xd positions);
private:
    const char* filename;
    int natoms;
    float timestep;
    float precision;
    matrix box;
    XDRFILE* xtc;
};

#endif // REPORTER_H