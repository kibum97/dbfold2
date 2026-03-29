#!/bin/bash
set -ex

# sm-run_namd_min.sh 
# Minimization script. Compatible with snakemake

echo sm-run_namd_min.sh; echo

if [ $# -ne 2 ]; then
    echo './sm-run_namd_min.sh inpdb output_directory'
    exit 1
fi

# module purge
# # module load centos6/0.0.1-fasrc01  ncf/1.0.0-fasrc01
# # module load hpc/namd-2.9
# module add GCC/7.3.0-2.30 OpenMPI/3.1.1
# module add NAMD/2.13-mpi


pdb=$1
output_directory=$2
: ${template:=/n/home00/vzhao/code/dbfold/MCPU/user_scripts/namd_minimization.template.conf}
: ${NAMDPATH:="/n/home00/vzhao/opt/NAMD_2.13_Linux-x86_64-multicore"}

fileroot=$(basename ${pdb} .pdb)
input_directory=$(dirname ${pdb})
coordinates=$(readlink -f "${input_directory}/${fileroot}.pdb")
structure=$(readlink -f "${input_directory}/${fileroot}.psf")

sed \
    -e s:XXXINPSF:"${structure}": \
    -e s:XXXINPDB:"${coordinates}": \
    -e s:XXXOUTROOT:"${fileroot}": \
    < "$template" \
    > "${output_directory}/namd_minimization_${fileroot}.conf"

${NAMDPATH}/namd2 ${output_directory}/namd_minimization_${fileroot}.conf
