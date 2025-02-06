#ifndef SYSTEM_H
#define SYSTEM_H

#include "topology.h"
#include "mcpu_forcefield.h"
#include "amino_acids.h"
#include "rotamer.h"

/**
 * @class System
 * @brief Represents a physical system with topology and force field.
 *
 * The System class encapsulates the topology and force field of a physical system.
 * It provides methods to initialize the system and paramters that are invariant thoroughout the simulation.
 * The parameters defined heres is mainly force field term that is parameterized based on the topology.
 */

struct TorsionIDData {
    size_t torsion_id;
    std::array<size_t, 4> torsion_atomIDs;
    std::vector<size_t> rotating_atomIDs;
};

struct AromaticIDData {
    size_t aromatic_id;
    std::array<size_t, 3> aromatic_atomIDs;
};

struct HBondIDData {
    size_t hbond_id;
    std::array<size_t, 7> donor_atomIDs;
    std::array<size_t, 8> acceptor_atomIDs;
};

class System {
public:
    System(Topology& topology, MCPUForceField& forcefield);
    ~System();

    // Refine topology by removing atoms
    // TODO: maybe split backbone and side chain positions?
    std::tuple<std::vector<size_t>, std::vector<size_t>> splitBackboneSidechain(Topology& topology);
    std::unordered_map<size_t, std::vector<size_t>> getRotatingAtomsMap(Topology& topology);
    std::unordered_map<size_t, std::vector<TorsionIDData>> getSideChainRotatingAtomsMap(Topology& topology, std::unordered_map<std::string, std::vector<TorsionData>> torsion_map);
    std::vector<AromaticIDData> getAromaticAtomsMap(Topology& topology, std::unordered_map<std::string, AromaticData> aromatic_map);
    std::vector<HBondIDData> getHBondAtomsMap(Topology& topology, HBondData hbond_map);
    std::unordered_map<size_t, std::vector<RotamerData>> getRotamerMap(Topology& topology, std::string rotamer_data_file);
    // Contact Energy - Mu Potential
    Eigen::MatrixXd getMuParameters() const;
    std::vector<int> getSmogType() const;
    // Hydrogen Bonding Energy
    // Torsional energy
    std::vector<BBTorsionParamArray> getBBTorsionParameters() const;
    std::vector<SCTorsionParamArray> getSCTorsionParameters() const;
    
private:
    Topology& topology; // Reference to Topology object
    MCPUForceField& forcefield; // Reference to ForceField object
    Eigen::MatrixXd muParamMatrix; // Matrix to store mu potential parameters
    std::vector<BBTorsionParamArray> bbTorsionParamVector; // Vector to store backbone torsion parameters
    std::vector<SCTorsionParamArray> scTorsionParamVector; // Vector to store side chain torsion parameters
    std::vector<int> smogTypeVector; // Vector to store Smog type
    int numAtoms; // Member variable to store the number of atoms
    int numResidues; // Member variable to store the number of residues

    void mupotential_parameters();
    void hbond_parameters();
    void torsion_parameters();

    int lookupAminoAcidIndex(const std::string& residueName);
};

#endif // SYSTEM_H
