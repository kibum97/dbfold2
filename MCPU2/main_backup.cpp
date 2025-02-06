#define FORCEFIELD 'MCPU'

#include <stdio.h>
#include <stdlib.h>
#include "mpi_util.h"
#include "read_config.h"
#include "topology.h"
#include "mcpu_forcefield.h"
#include "system.h"
#include "context.h"
#include "logger.h"
#include "load_pdb.h"


int main(int argc, char *argv[]) {
    int rank, size;

    // Initialize MPI
    mpi_initialize(&argc, &argv, &rank, &size);

    // Initialize Logger
    if (rank == 0) {
        Logger::initialize("simulation.log");
        Logger::log(INFO, "Simulation started");
    }

    // Check if configuration file is provided
    if (argc < 2) {
        if (rank == 0) {
            Logger::logf(ERROR, "Usage: %s <config_file>", argv[0]);
        }
        mpi_finalize();
        return EXIT_FAILURE;
    }

    // Read configuration file
    if (rank == 0) {
        SetProgramOptions(argv[1]);
        LogConfigurationDetails();
        
        // Initialize Topology and Forcefield
        Topology topology;
        topology.initialize(native_file);
        MCPUForceField forcefield;
        forcefield.initialize();

        std::cout << "Number of atoms: " << topology.getNumAtoms() << std::endl;

        // Initialize System
        System system(topology, forcefield);
        system.initialize();
        Eigen::MatrixXd muPotential = system.getMuParameters();

        // Initialize Context
        Context context(topology);
        context.initialize();
        Eigen::MatrixXd contacts1 = context.getContacts();
        std::cout << contacts1 << std::endl;

        // Compute contacts
        context.computeContacts();
        Eigen::MatrixXd contacts2 = context.getContacts();
        std::cout << contacts2 << std::endl;
        Eigen::MatrixXd binaryContacts = context.getBinaryContacts(5.0);
        std::cout << binaryContacts << std::endl;
        
        // Compute mu potential
        Eigen::MatrixXd result = muPotential.array() * binaryContacts.array();
        double mu_energy = result.sum();

        std::cout << "Sum of the result: " << mu_energy << std::endl;
    }

    // Initialize matrix

    // Finalize MPI
    mpi_finalize();

    // Finalize Logger
    if (rank == 0) {
        Logger::log(INFO, "Simulation ended");
        Logger::finalize();
    }

    return EXIT_SUCCESS;
}




/*
UNUSED CODE BLOCKS BACKUP

double sign(double a, double b) {
    if (b >= 0.)
        return fabs(a);
    else
        return -fabs(a);
}

double Context::computeDihedralAngle2(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c, const Eigen::Vector3d& d) {
    Eigen::Vector3d b1 = a - b;
    Eigen::Vector3d b2 = c - b;
    Eigen::Vector3d b3 = c - d;

    Eigen::Vector3d p = b1.cross(b2);
    Eigen::Vector3d q = b2.cross(b3);
    Eigen::Vector3d s = b3.cross(b1);

    
    double arg = p.dot(q) / sqrt(p.dot(p) * q.dot(q));
    arg = sign(std::min(fabs(arg), 1.0e0), arg);
    double angle = sign(acos(arg),s.dot(b2));

    return angle;
}
---------------------------------------

*/