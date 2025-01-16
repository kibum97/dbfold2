#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p shakhnovich,shared
#SBATCH -t 3-00:00
#SBATCH --mem=85000
#SBATCH -J test_min_%j                                                            
#SBATCH -o test_min_%j.out                                                         
#SBATCH -e test_min_%j.err                                                         
#SBATCH --mail-type=ALL                                                         
#SBATCH --mail-user=kibum_park@g.harvard.edu
# This file can be anywhere, but the dir that contains the cfg file
# must also have fold_potential_mpi. Then, the config_files dir should
# be one level up from that one.

module purge
module load gcc/12.2.0-fasrc01 openmpi/4.1.4-fasrc01

cd ./MCPU_min/1uao/MCPU_run/
srun -n ${SLURM_NTASKS} --mpi=pmi2 ./fold_potential_mpi ./cfg
