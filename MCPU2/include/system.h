#ifndef SYSTEM_H
#define SYSTEM_H

#include "topology.h"
#include "mcpu_forcefield.h"

/**
 * @class System
 * @brief Represents a physical system with topology and force field.
 *
 * The System class encapsulates the topology and force field of a physical system.
 * It provides methods to initialize the system and paramters that are invariant thoroughout the simulation.
 * The parameters defined heres is mainly force field term that is parameterized based on the topology.
 */

class Topology;
class MCPUForceField;

class System {
public:
    System(Topology& topology, MCPUForceField& forcefield);
    ~System();

    void initialize();
    // Contact Energy - Mu Potential
    Eigen::MatrixXd getMuParameters() const;
    std::vector<int> getSmogType() const;
    // Hydrogen Bonding Energy
    


private:
    Topology& topology; // Reference to Topology object
    MCPUForceField& forcefield; // Reference to ForceField object
    Eigen::MatrixXd muParamMatrix; // Matrix to store mu potential parameters
    std::vector<int> smogTypeVector; // Vector to store Smog type
    int numAtoms; // Member variable to store the number of atoms

    void mupotential_parameters();
    void hbond_parameters();
};

#endif // SYSTEM_H
