#!/bin/bash
set -ex

# run_psfgen.sh

# generates topology files for pdb files in a directory.


if [ $# -ne 2 ]; then
    echo './run_psfgen.sh structure_input_directory output_directory'
    exit 1
fi

module add vmd/1.9.3-fasrc01

input_directory=$1
output_directory=$2
: ${template:="/n/home00/vzhao/code/dbfold/MCPU/user_scripts/generate_topo.template.tcl"}

if [ ! -d "${output_directory}" ]; then
    mkdir -p "${output_directory}"
fi

for pdb in ${input_directory}/*pdb; do
    echo "Working on ${pdb}"
    fileroot=$(basename ${pdb} .pdb)
    
    sed \
	-e s:INPDBXXX:${pdb}: \
	-e s:OUTPSFXXX:"${output_directory}/${fileroot}.psf": \
	-e s:OUTPDBXXX:"${output_directory}/${fileroot}.pdb": \
	< "$template" \
	> generate_topo.tcl

    vmd -dispdev text -e generate_topo.tcl
done
