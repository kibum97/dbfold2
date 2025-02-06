#ifndef AMINO_ACIDS_H
#define AMINO_ACIDS_H

#include <unordered_map>
#include <vector>
#include <string>

inline const std::unordered_map<std::string, std::vector<std::pair<std::string, std::string>>> bondsTemplate = {
    {"ALA", {{"N", "CA"}, {"CA", "C"}, {"C", "O"}, {"CA", "CB"}}},
    {"GLY", {{"N", "CA"}, {"CA", "C"}, {"C", "O"}}},
    {"SER", {{"N", "CA"}, {"CA", "C"}, {"C", "O"}, {"CA", "CB"}, {"CB", "OG"}}},
    {"VAL", {{"N", "CA"}, {"CA", "C"}, {"C", "O"}, {"CA", "CB"}, {"CB", "CG1"}, {"CB", "CG2"}}},
    {"LEU", {{"N", "CA"}, {"CA", "C"}, {"C", "O"}, {"CA", "CB"}, {"CB", "CG"}, {"CG", "CD1"}, {"CG", "CD2"}}},
    {"ILE", {{"N", "CA"}, {"CA", "C"}, {"C", "O"}, {"CA", "CB"}, {"CB", "CG1"}, {"CB", "CG2"}, {"CG1", "CD1"}}},
    {"THR", {{"N", "CA"}, {"CA", "C"}, {"C", "O"}, {"CA", "CB"}, {"CB", "OG1"}, {"CB", "CG2"}}},
    {"CYS", {{"N", "CA"}, {"CA", "C"}, {"C", "O"}, {"CA", "CB"}, {"CB", "SG"}}},
    {"MET", {{"N", "CA"}, {"CA", "C"}, {"C", "O"}, {"CA", "CB"}, {"CB", "CG"}, {"CG", "SD"}, {"SD", "CE"}}},
    {"ASP", {{"N", "CA"}, {"CA", "C"}, {"C", "O"}, {"CA", "CB"}, {"CB", "CG"}, {"CG", "OD1"}, {"CG", "OD2"}}},
    {"ASN", {{"N", "CA"}, {"CA", "C"}, {"C", "O"}, {"CA", "CB"}, {"CB", "CG"}, {"CG", "OD1"}, {"CG", "ND2"}}},
    {"GLU", {{"N", "CA"}, {"CA", "C"}, {"C", "O"}, {"CA", "CB"}, {"CB", "CG"}, {"CG", "CD"}, {"CD", "OE1"}, {"CD", "OE2"}}},
    {"GLN", {{"N", "CA"}, {"CA", "C"}, {"C", "O"}, {"CA", "CB"}, {"CB", "CG"}, {"CG", "CD"}, {"CD", "OE1"}, {"CD", "NE2"}}},
    {"LYS", {{"N", "CA"}, {"CA", "C"}, {"C", "O"}, {"CA", "CB"}, {"CB", "CG"}, {"CG", "CD"}, {"CD", "CE"}, {"CE", "NZ"}}},
    {"ARG", {{"N", "CA"}, {"CA", "C"}, {"C", "O"}, {"CA", "CB"}, {"CB", "CG"}, {"CG", "CD"}, {"CD", "NE"}, {"NE", "CZ"}, {"CZ", "NH1"}, {"CZ", "NH2"}}},
    {"HIS", {{"N", "CA"}, {"CA", "C"}, {"C", "O"}, {"CA", "CB"}, {"CB", "CG"}, {"CG", "ND1"}, {"ND1", "CE1"}, {"CE1", "NE2"}, {"NE2", "CD2"}, {"CD2", "CG"}}},
    {"PHE", {{"N", "CA"}, {"CA", "C"}, {"C", "O"}, {"CA", "CB"}, {"CB", "CG"}, {"CG", "CD1"}, {"CD1", "CE1"}, {"CE1", "CZ"}, {"CZ", "CE2"}, {"CE2", "CD2"}, {"CD2", "CG"}}},
    {"TYR", {{"N", "CA"}, {"CA", "C"}, {"C", "O"}, {"CA", "CB"}, {"CB", "CG"}, {"CG", "CD1"}, {"CD1", "CE1"}, {"CE1", "CZ"}, {"CZ", "OH"}, {"CZ", "CE2"}, {"CE2", "CD2"}, {"CD2", "CG"}}},
    {"TRP", {{"N", "CA"}, {"CA", "C"}, {"C", "O"}, {"CA", "CB"}, {"CB", "CG"}, {"CG", "CD1"}, {"CD1", "NE1"}, {"NE1", "CE2"}, {"CE2", "CD2"}, {"CD2", "CG"}, {"CE2", "CZ2"}, {"CZ2", "CH2"}, {"CH2", "CZ3"}, {"CZ3", "CE3"}, {"CE3", "CD2"}}},
    {"PRO", {{"N", "CA"}, {"CA", "C"}, {"C", "O"}, {"CA", "CB"}, {"CB", "CG"}, {"CG", "CD"}, {"CD", "N"}}}
};

struct TorsionData {
    size_t torsion_id;
    std::vector<std::string> torsion_atoms;
    int num_rotating_atoms;
    std::vector<std::string> rotating_atoms;
};

inline const std::unordered_map<std::string, std::vector<TorsionData>> torsion_map = {
    {"ALA", {}},
    {"ARG", {{0, {"N", "CA", "CB", "CG"}, 6, {"CG", "CD", "NE", "CZ", "NH1", "NH2"}},
              {1, {"CA", "CB", "CG", "CD"}, 5, {"CD", "NE", "CZ", "NH1", "NH2"}},
              {2, {"CB", "CG", "CD", "NE"}, 4, {"NE", "CZ", "NH1", "NH2"}},
              {3, {"CG", "CD", "NE", "CZ"}, 3, {"CZ", "NH1", "NH2"}}}},
    {"ASN", {{0, {"N", "CA", "CB", "CG"}, 3, {"CG", "OD1", "ND2"}},
              {1, {"CA", "CB", "CG", "OD1"}, 2, {"OD1", "ND2"}}}},
    {"ASP", {{0, {"N", "CA", "CB", "CG"}, 3, {"CG", "OD1", "OD2"}},
              {1, {"CA", "CB", "CG", "OD1"}, 2, {"OD1", "OD2"}}}},
    {"CYS", {{0, {"N", "CA", "CB", "SG"}, 1, {"SG"}}}},
    {"GLN", {{0, {"N", "CA", "CB", "CG"}, 4, {"CG", "CD", "OE1", "NE2"}},
              {1, {"CA", "CB", "CG", "CD"}, 3, {"CD", "OE1", "NE2"}},
              {2, {"CB", "CG", "CD", "OE1"}, 2, {"OE1", "NE2"}}}},
    {"GLU", {{0, {"N", "CA", "CB", "CG"}, 4, {"CG", "CD", "OE1", "OE2"}},
              {1, {"CA", "CB", "CG", "CD"}, 3, {"CD", "OE1", "OE2"}},
              {2, {"CB", "CG", "CD", "OE1"}, 2, {"OE1", "OE2"}}}},
    {"GLY", {}},
    {"HIS", {{0, {"N", "CA", "CB", "CG"}, 5, {"CG", "ND1", "CD2", "CE1", "NE2"}},
              {1, {"CA", "CB", "CG", "ND1"}, 4, {"ND1", "CD2", "CE1", "NE2"}}}},
    {"ILE", {{0, {"N", "CA", "CB", "CG1"}, 3, {"CG1", "CG2", "CD1"}},
              {1, {"CA", "CB", "CG1", "CD1"}, 1, {"CD1"}}}},
    {"LEU", {{0, {"N", "CA", "CB", "CG"}, 3, {"CG", "CD1", "CD2"}},
              {1, {"CA", "CB", "CG", "CD1"}, 2, {"CD1", "CD2"}}}},
    {"LYS", {{0, {"N", "CA", "CB", "CG"}, 4, {"CG", "CD", "CE", "NZ"}},
              {1, {"CA", "CB", "CG", "CD"}, 3, {"CD", "CE", "NZ"}},
              {2, {"CB", "CG", "CD", "CE"}, 2, {"CE", "NZ"}},
              {3, {"CG", "CD", "CE", "NZ"}, 1, {"NZ"}}}},
    {"MET", {{0, {"N", "CA", "CB", "CG"}, 3, {"CG", "SD", "CE"}},
              {1, {"CA", "CB", "CG", "SD"}, 2, {"SD", "CE"}},
              {2, {"CB", "CG", "SD", "CE"}, 1, {"CE"}}}},
    {"PHE", {{0, {"N", "CA", "CB", "CG"}, 6, {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"}},
              {1, {"CA", "CB", "CG", "CD1"}, 5, {"CD1", "CD2", "CE1", "CE2", "CZ"}}}},
    {"PRO", {{0, {"N", "CA", "CB", "CG"}, 2, {"CG", "CD"}},
              {1, {"CA", "CB", "CG", "CD"}, 1, {"CD"}}}},
    {"SER", {{0, {"N", "CA", "CB", "OG"}, 1, {"OG"}}}},
    {"THR", {{0, {"N", "CA", "CB", "OG1"}, 2, {"OG1", "CG2"}}}},
    {"TRP", {{0, {"N", "CA", "CB", "CG"}, 9, {"CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"}},
              {1, {"CA", "CB", "CG", "CD1"}, 8, {"CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"}}}},
    {"TYR", {{0, {"N", "CA", "CB", "CG"}, 7, {"CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"}},
              {1, {"CA", "CB", "CG", "CD2"}, 6, {"CD1", "CD2", "CE1", "CE2", "CZ", "OH"}}}},
    {"VAL", {{0, {"N", "CA", "CB", "CG1"}, 2, {"CG1", "CG2"}}}},
};

struct AromaticData {
    std::vector<std::string> aromatic_ring_atoms;
};

inline const std::unordered_map<std::string, AromaticData> aromatic_map = {
    {"PHE", {{"CG", "CE1", "CE2"}}},
    {"TRP", {{"CG", "CZ2", "CZ3"}}},
    {"RING", {{"CG", "CE1", "CE2"}}},
    //{"TYR", {{{"CG", "CE1", "CE2"}}}}
};

struct HBondData {
    std::vector<std::pair<int, std::string>> donor_atoms;
    std::vector<std::pair<int, std::string>> acceptor_atoms;
};

inline const HBondData hbond_map = {
    { // Donor
        {{0, "N"}, {0, "CA"}, {-1, "C"}, {0, "C"}, {-1, "CA"}, {-1, "N"}, {1, "CA"}}
    }, 
    { // Acceptor
        {{0, "O"}, {0, "C"}, {0, "CA"}, {1, "N"}, {-1, "N"}, {1, "CA"}, {1, "C"}, {-1, "CA"}}
    }
};

inline const std::unordered_map<std::string, size_t> residueType_map = {
    {"ALA", 0}, {"ARG", 1}, {"ASN", 2}, {"ASP", 3}, {"CYS", 4},
    {"GLN", 5}, {"GLU", 6}, {"GLY", 7}, {"HIS", 8}, {"ILE", 9},
    {"LEU", 10}, {"LYS", 11}, {"MET", 12}, {"PHE", 13}, {"PRO", 14},
    {"SER", 15}, {"THR", 16}, {"TRP", 17}, {"TYR", 18}, {"VAL", 19}
};

#endif // AMINO_ACIDS_H