#!/bin/bash
set -ex

# run_namd_min.sh 
# Minimization for multiple structures in one directory. single core

echo run_namd_min.sh; echo

if [ $# -ne 2 ]; then
    echo './run_namd_min.sh structure_input_directory output_directory'
    exit 1
fi

# module purge
# # module load centos6/0.0.1-fasrc01  ncf/1.0.0-fasrc01
# # module load hpc/namd-2.9
# module add GCC/7.3.0-2.30 OpenMPI/3.1.1
# module add NAMD/2.13-mpi


input_directory=$1
output_directory=$2
: ${template:=/n/home00/vzhao/opt/MCPU/user_scripts/namd_minimization.template.conf}
: ${NAMDPATH:="/n/home00/vzhao/opt/NAMD_2.13_Linux-x86_64-multicore"}

# SLURM associated
: ${PARTITION:="shakhnovich,shakgpu,shared"}
: ${TIME:=60}
: ${DRYRUN:=}

if [ ! -d "${output_directory}" ]; then
    mkdir -p "${output_directory}"
fi

for pdb in ${input_directory}/*pdb; do
    echo "Working on ${pdb}"
    fileroot=$(basename ${pdb} .pdb)
    coordinates=$(readlink -f "${input_directory}/${fileroot}.pdb")
    structure=$(readlink -f "${input_directory}/${fileroot}.psf")
    
    sed \
	-e s:XXXINPSF:"${structure}": \
	-e s:XXXINPDB:"${coordinates}": \
	-e s:XXXOUTROOT:"${fileroot}": \
	< "$template" \
	> "${output_directory}/namd_minimization_${fileroot}.conf"

    $([ -n "$DRYRUN" ] && echo echo) \
    	sbatch -p ${PARTITION} -n 1 -t ${TIME} --mem 512 \
    	-o ${output_directory}/${fileroot}.slurm -J "min_${fileroot}" \
    	--wrap "\"${NAMDPATH}/namd2\" ${output_directory}/namd_minimization_${fileroot}.conf"
    # "${NAMDPATH}/namd2" ${output_directory}/namd_minimization_${fileroot}.conf
    sleep 1
done
