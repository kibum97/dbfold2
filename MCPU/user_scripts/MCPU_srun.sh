#!/bin/bash
#SBATCH --no-requeue
# This file can be anywhere, but the dir that contains the cfg file
# must also have fold_potential_mpi. Then, the config_files dir should
# be one level up from that one.

module purge
module load gcc/9.3.0-fasrc01 openmpi/4.0.2-fasrc01

# srun -n ${SLURM_NTASKS} --mpi=pmix fold_potential_mpi "$@"
# orterun -n ${SLURM_NTASKS} fold_potential_mpi "$@"

gcc -c /n/home01/kibumpark/pkg/dbfold/MCPU/src_mpi_umbrella/rng.c -o /n/home01/kibumpark/pkg/dbfold/MCPU/src_mpi_umbrella/rng.o
mpicc -O3 -o /n/home01/kibumpark/pkg/dbfold/MCPU/src_mpi_umbrella/fold_potential_mpi /n/home01/kibumpark/pkg/dbfold/MCPU/src_mpi_umbrella/backbone.c -lm /n/home01/kibumpark/pkg/dbfold/MCPU/src_mpi_umbrella/rng.o

cfg="$(readlink -f "$1")"
cd "$(dirname "$cfg")"

srun -n ${SLURM_NTASKS} --mpi=pmix ./fold_potential_mpi "$(basename "$cfg")"
# else
#     srun -n ${SLURM_NTASKS} --mpi=pmi2 ./fold_potential_mpi "$(basename "$cfg")"
# fi
