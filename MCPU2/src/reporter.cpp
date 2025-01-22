#include "reporter.h"

Reporter::Reporter(const char* filename, int natoms, float timestep, float precision, const Eigen::Matrix3d& periodicBox)
    : filename(filename), natoms(natoms), timestep(timestep), precision(precision) {
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            this->box[i][j] = static_cast<float>(periodicBox(i, j));
        }
    }
    xtc = xdrfile_open(filename, "w");
    if (!xtc) {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
    }
}

Reporter::~Reporter() {
    if (xtc) {
        xdrfile_close(xtc);
    }
}

bool Reporter::writeTrajectory(int frame_idx, Eigen::Matrix3Xd positions) {
    if (!xtc) return false;
    rvec* x = new rvec[natoms];
    for (int atom_idx = 0; atom_idx < natoms; ++atom_idx) {
        x[atom_idx][0] = positions(0, atom_idx);
        x[atom_idx][1] = positions(1, atom_idx);
        x[atom_idx][2] = positions(2, atom_idx);
    }

    if (write_xtc(xtc, natoms, frame_idx, timestep * frame_idx, box, x, precision) != exdrOK) {
        std::cerr << "Error: Could not write frame " << frame_idx << " to file." << std::endl;
        delete[] x;
        return false;
    }

    delete[] x;
 
    return true;
}