#!/bin/bash
set -ex

<<'EOF'

After NAMD 5000 step minimization on the structures,
1. Create from template a TCL script to extract minimized structure
2. Call vmd to run script
3. Take output temp pdb and run it through sed to change some aliased names

THIS SCRIPT IS FRAGILE. ONLY ONE INSTANCE CAN BE RUN AT A TIME.

EOF

echo extract_minimized_structure.sh; echo

if [ $# -ne 3 ]; then
    echo "./extract_minimized_structure.sh structure_dir trajectory_dir output_dir"
    exit 1
fi

module purge
module add vmd/1.9.3-fasrc01

structure_directory=$1
trajectory_directory=$2
output_directory=$3
template=/n/home00/vzhao/code/dbfold/MCPU/user_scripts/extract_pdb0.tcl.template

if [ ! -d "${output_directory}" ]; then
    mkdir -p "${output_directory}"
fi

for pdb in ${structure_directory}/*pdb; do
    rm -rf MIN_OUT_PDB_TEMP.pdb extract_minimized_structure.tcl
    fileroot=$(basename ${pdb} .pdb)
    echo "Working on ${fileroot}"
    # coordinates=${structure_directory}/${fileroot}.pdb
    structure=${structure_directory}/${fileroot}.psf
    trajectory=${trajectory_directory}/${fileroot}.dcd
    
    sed \
	-e s:INPSFXXX:"${structure}": \
	-e s:INDCDXXX:"${trajectory}": \
	< "$template" \
	> extract_minimized_structure.tcl
    if [ ! -f "extract_minimized_structure.tcl" ]; then
	echo "no extract_minimized_structure.tcl"
	echo "sed failed?"
	exit 1
    fi

    vmd -dispdev text -e extract_minimized_structure.tcl

    if [ -f MIN_OUT_PDB_TEMP.pdb ]; then
	sed \
	    -e s/HSE/HIS/ \
	    -e "/ILE/ s/CD /CD1/" \
	    < MIN_OUT_PDB_TEMP.pdb \
	    > ${output_directory}/${fileroot}_minimized.pdb
    else
	echo "no MIN_OUT_PDB_TEMP.pdb"
	echo "extract must have failed"
	exit 1
    fi
done
