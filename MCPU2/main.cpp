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
#include "amino_acids.h"
#include "rotamer.h"
#include <chrono>
#include <random>

/*
TODO
- implement hydrogen bond energy and aromatic energy and torsional energy
- implement the move acceptance criteria
- implement the MC steps
- implement the replica exchange
- Make reading files part in the same style
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
        std::unordered_map<size_t, std::vector<size_t>> rotatingBBAtomsMap = system.getRotatingAtomsMap(topology);
        std::unordered_map<size_t, std::vector<TorsionIDData>> rotatingSCAtomsMap = system.getSideChainRotatingAtomsMap(topology, torsion_map);
        std::vector<AromaticIDData> aromaticAtomsMap = system.getAromaticAtomsMap(topology, aromatic_map);
        for (const auto& [key, value] : rotatingSCAtomsMap) {
            std::cout << key << ": ";
            for (const auto& data : value) {
                std::cout << data.torsion_atomIDs[0] << " " << data.torsion_atomIDs[1] << " " << data.torsion_atomIDs[2] << " " << data.torsion_atomIDs[3] << " ";
            }
            std::cout << std::endl;
        }
        std::vector<HBondIDData> hbondAtomsMap = system.getHBondAtomsMap(topology, hbond_map);
        std::cout << "Aromatic atoms map: " << std::endl;
        std::cout << aromaticAtomsMap.size() << std::endl;
        std::cout << "HBond atoms map: " << hbondAtomsMap.size() << std::endl;
        std::unordered_map<size_t, std::vector<RotamerData>> rotamerIDMap = system.getRotamerMap(topology, rotamer_data_file);
        std::cout << "Rotating atoms map: " << std::endl;
        for (const auto& [key, value] : rotatingBBAtomsMap) {
            std::cout << key << ": ";
            for (const auto& atomID : value) {
                std::cout << atomID << " ";
            }
            std::cout << std::endl;
        }
        std::cout << "System initialized successfully." << std::endl;
        Eigen::MatrixXd muPotential = system.getMuParameters();
        std::vector<BBTorsionParamArray> bbTorsionParamVector = system.getBBTorsionParameters();
        std::vector<SCTorsionParamArray> scTorsionParamVector = system.getSCTorsionParameters();
        std::map<int, double> hbond_map = forcefield.getHBondPotential();
        std::map<int, double> aromatic_map = forcefield.getAromaticPotential();
        HBondParamArray seq_based_hbondPotential = forcefield.getSequenceDependentHBondPotential();

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
        Context context(topology, positions, hbondAtomsMap, aromaticAtomsMap);
        std::cout << "Context initialized successfully." << std::endl;
        context.computeSidechainAngles(rotatingSCAtomsMap);
        std::vector<std::array<double, 4>> sidechain_torsions = context.getSideChainTorsion();
        std::cout << "Sidechain torsion initialized" << std::endl;
        context.computeHBondStates();
        std::cout << "HBond states initialized" << std::endl;
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
        context.initializeContext();
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
        Eigen::Vector3d axis;
        std::cout << "Size of bb_atomids: " << bb_atomids.size() << std::endl;
        //std::exit(1);
        for (int i = 0; i < bb_atomids.size()-1; ++i) {
            // Create a deep copy of the context object
            Context new_context = context;
            auto start = std::chrono::high_resolution_clock::now(); // Start time

            //std::cout << "Performing concerted rotation..." << "Step " << i << std::endl;
            //std::cout << "Performing partial rotation..." << std::endl;
            //std::cout << "Atom 1: " << bb_atomids[i] << " Atom 2: " << bb_atomids[i+1] << std::endl;
            //std::cout << "Atom 1: " << topology.atoms[bb_atomids[i]].atomName << " Atom 2: " << topology.atoms[bb_atomids[i+1]].atomName << std::endl;
            axis = positions.col(bb_atomids[i]) - positions.col(bb_atomids[i+1]);
            axis.normalize();
            //std::cout << "Axis vector: " << axis << std::endl;
            movedAtomIDs = rotatingBBAtomsMap[bb_atomids[i]];
            std::vector<size_t> movedResidueIDs;
            size_t startResID = topology.atoms[bb_atomids[i]].resID;
            //std::cout << "Start ResID: " << startResID << "/" << topology.getNumResidues() <<std::endl;
            for (size_t resID = startResID; resID < topology.getNumResidues(); ++resID) {
                movedResidueIDs.push_back(resID);
            }
            std::cout << "Rotation will happen" << std::endl;
            positions = integrator.partialRotation(positions, movedAtomIDs, 90.0, axis);
            std::cout << "BB Rotation done" << std::endl;
            positions = integrator.sidechainRotation(startResID, positions, rotamerIDMap, rotatingSCAtomsMap, sidechain_torsions);
            std::cout << "SC Rotation done" << std::endl;
            //std::cout << "Number of moved residues: " << movedResidueIDs.size() << std::endl;
            auto part1_end = std::chrono::high_resolution_clock::now();
            new_context.updateContext(movedAtomIDs, movedResidueIDs, positions);
            //std::cout << "Context updated" << std::endl;
            //reporter.writeTrajectory(i+1, positions);
            // Compute New energy
            //std::cout << "Computing contacts..." << std::endl;
            auto part2_end = std::chrono::high_resolution_clock::now();
            new_context.computeEnergy(muPotential, bbTorsionParamVector, scTorsionParamVector, hbond_map, seq_based_hbondPotential, aromatic_map);
            std::cout << "Energy computed" << std::endl;
            auto part3_end = std::chrono::high_resolution_clock::now();
            // Perform MH step
            // Perform Metropolis-Hastings step
            if (new_context.getEnergy().total_energy < context.getEnergy().total_energy) {
                std::swap(context, new_context);
                std::cout << "New state accepted" << std::endl;
            } else {
                double acceptance_probability = exp(-(new_context.getEnergy().total_energy - context.getEnergy().total_energy) / 0.4);
                std::random_device rd;
                std::mt19937 gen(rd());
                std::uniform_real_distribution<> dis(0.0, 1.0);
                double randomValue = dis(gen);
                if (randomValue < acceptance_probability) {
                    std::swap(context, new_context);
                    std::cout << "New state accepted" << std::endl;
                } else {
                    std::cout << "New state rejected" << std::endl;
                }
            }
            auto part4_end = std::chrono::high_resolution_clock::now();
            
            // Calculate elapsed time for each part
            std::chrono::duration<double> part1_elapsed = part1_end - start;
            std::chrono::duration<double> part2_elapsed = part2_end - part1_end;
            std::chrono::duration<double> part3_elapsed = part3_end - part2_end;
            std::chrono::duration<double> part4_elapsed = part4_end - part3_end;
            
            // Print elapsed time for each part
            std::cout << "Elapsed time for step " << i << " - Part 1: " << part1_elapsed.count() << " seconds" << std::endl;
            std::cout << "Elapsed time for step " << i << " - Part 2: " << part2_elapsed.count() << " seconds" << std::endl;
            std::cout << "Elapsed time for step " << i << " - Part 3: " << part3_elapsed.count() << " seconds" << std::endl;
            std::cout << "Elapsed time for step " << i << " - Part 4: " << part4_elapsed.count() << " seconds" << std::endl;
            
            auto end = std::chrono::high_resolution_clock::now(); // End time
            std::chrono::duration<double> total_elapsed = end - start; // Calculate total elapsed time
            std::cout << "Total elapsed time for step " << i << ": " << total_elapsed.count() << " seconds" << std::endl;
        }
        std::cout << "Number of residues: " << topology.getNumResidues() << std::endl;

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

    // Finalize MPI
    mpi_finalize();

    // Finalize Logger
    if (rank == 0) {
        Logger::log(INFO, "Simulation ended");
        Logger::finalize();
    }

    return EXIT_SUCCESS;
}