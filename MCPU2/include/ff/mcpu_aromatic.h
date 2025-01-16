#ifndef MCPU_AROMATIC_H
#define MCPU_AROMATIC_H

#include "mcpu_forcefield.h"
#include <unordered_map>
#include <vector>
#include <string>

std::unordered_map<int, int> InitAromatic(const std::string& filename);
std::vector<int> AromaticsFromTopology(const std::vector<Atom>& atoms);

#endif // MCPU_AROMATIC_H