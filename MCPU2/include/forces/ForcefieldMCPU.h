#ifndef FORCEFIELD_MCPU_H
#define FORCEFIELD_MCPU_H

#define NUM_SMOG_TYPES 84

#include <Eigen/Dense>
#include <map>
#include <vector>

class ProteinStructure;

using BackboneTorsionPerTriplet = std::array<std::array<std::array<std::array<int, 6>, 6>, 6>, 6>;
using BackboneTorsionParams =
    std::array<std::array<std::array<BackboneTorsionPerTriplet, 20>, 20>, 20>;
using SidechainTorsionPerTriplet =
    std::array<std::array<std::array<std::array<int, 12>, 12>, 12>, 12>;
using SidechainTorsionParams =
    std::array<std::array<std::array<SidechainTorsionPerTriplet, 20>, 20>, 20>;
using HBondParams          = std::array<std::array<std::array<double, 20>, 20>, 3>;
using ContactCutoff2Matrix = Eigen::Matrix<float, NUM_SMOG_TYPES, NUM_SMOG_TYPES>;

class ForcefieldMCPU {
   public:
    ForcefieldMCPU();
    ~ForcefieldMCPU();

    // Member variables
    Eigen::MatrixXf                                        muPotential;
    BackboneTorsionParams                                  bbTorsions;
    SidechainTorsionParams                                 scTorsions;
    std::map<int, float>                                   hbondPotential;
    HBondParams                                            sequenceDependentHBondPotential;
    std::map<int, float>                                   aromaticPotential;
    std::unordered_map<std::string, std::pair<int, float>> smogTypeMap;
    ContactCutoff2Matrix                                   contact_cutoff2;

   private:
    // Smog Type
    // std::tuple<std::vector<int>, std::vector<double>, std::vector<int>>
    // assignSmogType(ProteinStructure& protein);
    // Contact Energy - Mu Potential
    void loadMuPotential(const std::string &filename);
    void loadSmogTypeMap(const std::string &filename);
    void initialize_contact_cutoff2();
    // Hydrogen Bonding Energy
    void loadHbondPotential(const std::string &filename);
    void loadSequenceDependentHBondPotential(const std::string &filename);
    // Torsion Energy - Backbone and Sidechain
    void loadTorsionalPotential(const std::string &triplet_file, const std::string &sctorsion_file);
    // Aromatic Energy
    void loadAromaticPotential(const std::string &filename);
};

#endif  // FORCEFIELD_MCPU_H