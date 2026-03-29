#!/bin/bash
module purge
module load gcc openmpi cmake

cd /n/home01/kibumpark/pkg/dbfold2/MCPU2
cmake -DCMAKE_CXX_FLAGS="-O3 -march=native -funroll-loops -flto" --preset=default
#cmake -DCMAKE_CXX_FLAGS="-O3" --preset=default
#cmake --preset=default
cmake --build build --verbose

#srun -n 1 --mpi=pmi2 ./build/MCPU2 /n/home01/kibumpark/pkg/dbfold2/examples/abce1/config.yaml
