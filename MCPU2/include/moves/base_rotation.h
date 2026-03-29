#ifndef BASE_ROTATION_H
#define BASE_ROTATION_H

#include <Eigen/Dense>
#include <map>
#include <vector>

class BaseRotation {
   public:
    virtual void rotate(int pivot_atom) const = 0;
    virtual void generate_rotating_atom_map() = 0;
    virtual ~BaseRotation()                   = default;

   protected:
    std::map<int, std::vector<int>>
        rotating_atom_map;  // Map of atom indices to their corresponding rotating atoms
};

#endif  // BASE_ROTATION_H