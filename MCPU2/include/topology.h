#ifndef TOPOLOGY_H
#define TOPOLOGY_H

#include <string>
#include <vector>
#include <Eigen/Dense>

class Atom
{
public:
    std::string atomName; // Atom name (e.g., "O", "CA")
    std::string element;  // Element (e.g., "O", "C")
    int atomNum;          // Atom number
    size_t atomID;        // Atom index in the vector of atoms
    size_t resID;         // Residue index in the vector of residues

    Atom(const std::string& atomName, const std::string& element, int atomNum, size_t atomID, size_t resID)
        : atomName(atomName), element(element), atomNum(atomNum), atomID(atomID), resID(resID) {}
};

class Residue
{
public:
    std::string resName;     // Residue name (e.g., "TYR")
    int resNum;              // Residue sequence number
    size_t resID;            // Residue index in the vector of residues
    std::vector<Atom> atoms; // Vector of atoms in residue
    size_t chainID;          // Chain index in the vector of chains

    Residue(const std::string& resName, int resNum, size_t resID, size_t chainID, const std::vector<Atom>& atoms)
        : resName(resName), resNum(resNum), resID(resID), chainID(chainID), atoms(atoms) {}
};

class Chain
{
public:
    char chainName;                // Chain ID (e.g., 'A')
    size_t chainID;                // Chain index in the vector of chains
    std::vector<Residue> residues; // Vector of residues in chain

    Chain(char chainName, size_t chainID, const std::vector<Residue>& residues)
        : chainName(chainName), chainID(chainID), residues(residues) {}
};

class Bond
{
public:
    Atom atom1;
    Atom atom2;
};

class Dihedral
{
public:
    Atom atom1;
    Atom atom2;
    Atom atom3;
    Atom atom4;
};

class Hbond
{
public:
    std::vector<Atom> donor_atoms;
    std::vector<Atom> acceptor_atoms;
};

class Aromatic
{
public:
    std::vector<Atom> atoms;
};

class Topology
{
public:
    int getNumAtoms() const;
    int getNumResidues() const;
    int getNumChains() const;

    void setPeriodicBox(const Eigen::Matrix3d& periodicBox);
    Eigen::Matrix3d getPeriodicBox() const;

    std::vector<Atom> atoms;       // Member variable to store atoms
    std::vector<Residue> residues; // Member variable to store residues; DO I NEED THIS?
    std::vector<Chain> chains;     // Member variable to store chains

    void removeAtomsByID(const std::vector<int>& atomIDs);
    void updateTopology();
    void printTopology();

private:
    Eigen::Matrix3d periodicBox; // Member variable to store box size
    int numAtoms;            // Member variable to store number of atoms
    int numResidues;         // Member variable to store number of residues
    int numChains;           // Member variable to store number of chains
};

#endif // TOPOLOGY_H