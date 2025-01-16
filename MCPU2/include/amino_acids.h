#include <unordered_map>
#include <vector>
#include <string>

const std::unordered_map<std::string, std::vector<std::pair<std::string, std::string>>> bondsTemplate = {
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

