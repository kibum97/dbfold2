#!/bin/bash
module purge
module load gcc openmpi cmake
#module load intel/25.0.1-fasrc01
#module load Mambaforge intelmpi cmake

cd /n/home01/kibumpark/pkg/dbfold2/MCPU2
#cmake -DCMAKE_CXX_FLAGS="-O3" --preset=default
cmake --preset=default
cmake --build build --verbose

#srun -n 1 --mpi=pmi2 ./build/MCPU2 ../examples/chignolin/config.yaml 
