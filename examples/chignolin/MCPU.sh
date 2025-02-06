#!/bin/bash
module purge
module load gcc openmpi cmake

cd /n/home01/kibumpark/pkg/dbfold2/MCPU2
cmake -DCMAKE_CXX_FLAGS="-O3" --preset=default
cmake --build build

#srun -n 1 --mpi=pmi2 ./build/MCPU2 ../examples/chignolin/config.yaml 
