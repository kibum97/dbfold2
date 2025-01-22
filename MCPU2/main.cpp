#define FORCEFIELD 'MCPU'

#include <stdio.h>
#include <stdlib.h>
#include "mpi_util.h"
#include <iomanip>
#include "read_config.h"
#include "topology.h"
#include "logger.h"
#include "load_pdb.h"
#include "context.h"
#include "mcpu_forcefield.h"
#include "system.h"
#include "integrator.h"
#include "reporter.h"

/*
TODO
- position updates should also applied to context (update context after every move?)
- define what moves which atoms
- implement hydrogen bond energy and aromatic energy and torsional energy
- implement the move acceptance criteria
- implement the MC steps
- implement the replica exchange
- implement the output
*/


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

    if (rank == 0) {
        // Read configuration file
        SetProgramOptions(argv[1]);
        LogConfigurationDetails();

        // Initialize Topology
        // At this stage topology will be simply read from a PDB file
        // without considering the forcefield
        std::cout << "Parsing PDB file..." << std::endl;
        auto [pdb_topology, pdb_positions] = parsePDB(native_file);
        std::cout << "Parsed PDB file successfully." << std::endl;

        // Initialize System
        // At this stage, forcefield will be considered
        // and topology will be refined based on the forcefield
        MCPUForceField forcefield;
        // Get Smog Type
        auto [smogTypeVector, remove_atom_ids] = forcefield.getSmogType(pdb_topology);
        // Refine topology
        std::cout << remove_atom_ids.size() << std::endl;
        auto [topology, positions] = forcefield.removeAtomsByID(pdb_topology, pdb_positions, remove_atom_ids); // TODO: Clean up the function
        std::cout << "Initializing system..." << std::endl;
        System system(topology, forcefield);
        auto [bb_atomids, sc_atomids] = system.splitBackboneSidechain(topology);
        std::cout << "System initialized successfully." << std::endl;
        Eigen::MatrixXd muPotential = system.getMuParameters();

        // Initialize reporter
        const char* filename = "trajectory.xtc";
        int natoms = topology.getNumAtoms();
        int nframes = 100;
        float timestep = 1;
        float precision = 1000.0f;
        Eigen::Matrix3d periodicBox = topology.getPeriodicBox();
        Reporter reporter(filename, natoms, timestep, precision, periodicBox);

        // Initialize Context
        // Context contains the values that changes over MC steps
        std::cout << "Initializing context..." << std::endl;
        Context context(topology, positions);
        std::cout << "Context initialized successfully." << std::endl;
        // Assign atom IDs to cells
        CellList cellList = context.getCellList();
        /*
        for (int i = 0; i < cellList.cells.size(); ++i) {
            std::cout << "Atom IDs in cell " << i << ": ";
            for (const auto& atomID : cellList.cells[i].atomIDs) {
                std::cout << atomID << " ";
            }
            std::cout << std::endl;
        }
        */
        
        reporter.writeTrajectory(0, positions);

        // Compute contacts
        Eigen::MatrixXd contacts = context.getContacts();
        std::cout << "Computing contacts..." << std::endl;
        for (int i = 0; i < cellList.cells.size(); ++i) {
            context.computeContacts(i);
        }
        contacts = context.getContacts();
        std::cout << "Contacts computed successfully." << std::endl;
        Eigen::MatrixXd binaryContacts = context.getBinaryContacts(5.0);
        // Compute Mu potential energy       
        Eigen::MatrixXd result = muPotential.array() * binaryContacts.array();
        double mu_energy = result.sum();
        std::cout << "Mu potential energy: " << mu_energy << std::endl;

        // Make movements
        Integrator integrator;
        std::vector<size_t> movedAtomIDs;
        for (size_t i = 4; i < topology.getNumAtoms(); ++i) {
            movedAtomIDs.push_back(i);
        }
        Eigen::Vector3d axis;
        std::cout << "Performing partial rotation..." << std::endl;
        std::cout << "Atom 1: " << bb_atomids[3] << " Atom 2: " << bb_atomids[4] << std::endl;
        std::cout << "Atom 1: " << topology.atoms[bb_atomids[3]].atomName << " Atom 2: " << topology.atoms[bb_atomids[4]].atomName << std::endl;
        axis = positions.col(bb_atomids[3]) - positions.col(bb_atomids[4]);
        axis.normalize();
        std::cout << "Axis vector: " << axis << std::endl;
        positions = integrator.partialRotation(positions, movedAtomIDs, 90.0, axis);
        context.updateContext(movedAtomIDs, positions);
        reporter.writeTrajectory(1, positions);
        
        // Compute New energy
        std::cout << "Computing contacts..." << std::endl;
        contacts = context.getContacts();
        std::cout << "Contacts computed successfully." << std::endl;
        binaryContacts = context.getBinaryContacts(5.0);
        // Compute Mu potential energy       
        result = muPotential.array() * binaryContacts.array();
        mu_energy = result.sum();
        std::cout << "Mu potential energy: " << mu_energy << std::endl;

        // Save to pdb
        std::ofstream output("output.pdb");
        for (int i = 0; i < topology.getNumAtoms(); ++i) {
            auto atom = topology.atoms[i];
            auto residue = topology.residues[atom.resID];
            auto chain = topology.chains[residue.chainID];
            output << "ATOM  " << std::setw(5) << i + 1 << " " << std::setw(4) << atom.atomName << " " << residue.resName << " " << chain.chainName << std::setw(4) << residue.resNum << "    ";
            output << std::fixed << std::setprecision(3) << std::setw(8) << positions(0, i) << std::setw(8) << positions(1, i) << std::setw(8) << positions(2, i) << std::endl;
        }
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