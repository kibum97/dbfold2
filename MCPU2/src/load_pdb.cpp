#include "load_pdb.h"

std::tuple<Topology, Eigen::Matrix3Xd> parsePDB(const std::string &filename) {
    Topology topology;
    std::vector<Eigen::Vector3d> all_coords;
    Eigen::Matrix3Xd positions;
    std::ifstream pdbFile(filename);
    if (!pdbFile.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        exit(1);
    }

    std::unordered_map<char, size_t> chainMap; // Map for chain ID to index
    std::unordered_map<int, size_t> residueMap; // Map for residue sequence to index within a chain
    std::string line;

    while (std::getline(pdbFile, line)) {
        // Process only ATOM or HETATM records
        if (line.substr(0, 6) == "ATOM  " || line.substr(0, 6) == "HETATM") {
            // Extract atom information based on fixed-column positions in the PDB format
            int atomNum = std::stoi(line.substr(6, 5)); // Atom serial number (columns 7-11)
            std::string atomName = trim(line.substr(12, 4));  // Atom name (columns 13-16)
            std::string element = trim(line.substr(76, 2));   // Element (columns 77-78)
            std::string resName = trim(line.substr(17, 3));   // Residue name (columns 18-20)
            char chainName = line[21];           // Chain ID (column 22)
            int resNum = std::stoi(line.substr(22, 4)); // Residue sequence number (columns 23-26)
            Eigen::Vector3d coord;
            coord(0) = std::stod(line.substr(30, 8)); // X coordinate (columns 31-38)
            coord(1) = std::stod(line.substr(38, 8));      // Y coordinate (columns 39-46)
            coord(2) = std::stod(line.substr(46, 8));      // Z coordinate (columns 47-54)
            
            // Check if chain ID is already in the map
            size_t chainIndex;
            if (chainMap.find(chainName) == chainMap.end()) {
                topology.chains.push_back(Chain{chainName, topology.chains.size(), {}});
                chainMap[chainName] = topology.chains.size() - 1;
            }
            chainIndex = chainMap[chainName];

            // Check if residue sequence number is already in the map
            size_t residueIndex;
            if (residueMap.find(resNum) == residueMap.end()) {
                topology.chains[chainIndex].residues.push_back(Residue{resName, resNum, topology.chains[chainIndex].residues.size(), chainIndex, {}});
                residueMap[resNum] = topology.chains[chainIndex].residues.size() - 1;
            }
            residueIndex = residueMap[resNum];

            // Add atom to residue
            topology.chains[chainIndex].residues[residueIndex].atoms.push_back(Atom{atomName, element, atomNum, topology.atoms.size(), residueIndex});

            // Add coordinates to all_coords
            all_coords.push_back(coord);
        }
    }

    // Merge coordinates into a single matrix and set periodic box
    positions = Eigen::Matrix3Xd::Zero(3, all_coords.size());
    for (int i = 0; i < all_coords.size(); ++i) {
        positions.col(i) = all_coords[i];
    }
    
    // Calculate the minimum and maximum coordinates along each axis
    Eigen::Vector3d minCoord = positions.rowwise().minCoeff();
    Eigen::Vector3d maxCoord = positions.rowwise().maxCoeff();

    // Define the periodic box vectors based on the min and max coordinates
    Eigen::Vector3d a = Eigen::Vector3d(maxCoord(0) - minCoord(0), 0.0, 0.0);
    Eigen::Vector3d b = Eigen::Vector3d(0.0, maxCoord(1) - minCoord(1), 0.0);
    Eigen::Vector3d c = Eigen::Vector3d(0.0, 0.0, maxCoord(2) - minCoord(2));

    // Set the periodic box in the topology
    Eigen::Matrix3d periodicBox;
    periodicBox.col(0) = a;
    periodicBox.col(1) = b; 
    periodicBox.col(2) = c;
    topology.setPeriodicBox(periodicBox);

    return std::make_tuple(topology, positions);
}

std::string trim(const std::string &str) {
    // Trim whitespace from a string
    size_t start = str.find_first_not_of(" \t");
    size_t end = str.find_last_not_of(" \t");
    return (start == std::string::npos) ? "" : str.substr(start, end - start + 1);
}