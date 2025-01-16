#ifndef LOAD_PDB_H
#define LOAD_PDB_H

#include "topology.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <Eigen/Dense>

std::tuple<Topology, Eigen::Matrix3Xd> parsePDB(const std::string &filename);
std::string trim(const std::string &str);

#endif // LOAD_PDB_H
