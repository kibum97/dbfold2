#ifndef CONCERTED_ROTATION_H
#define CONCERTED_ROTATION_H

#include <Eigen/Dense>
#include <unordered_map>
#include <vector>

class SystemMCPU;

class ConcertedRotation {
   public:
    ConcertedRotation(SystemMCPU &system);

    static Eigen::Matrix3Xf rotate(Eigen::Matrix3Xf &positions, int &chain_id, int &residue_id,
                                   float angle);
    void generate_pivot_atom_id_to_rotating_atom_ids_map(const std::vector<int> &backbone_atom_ids);

    static std::vector<std::unordered_map<int, std::vector<int>>> rotation_atom_ids;
    static int                                                    n_solutions;

   private:
    SystemMCPU &system_;  // Reference to the system
};

#endif  // CONCERTED_ROTATION_H