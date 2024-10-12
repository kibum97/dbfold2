#!/bin/bash
set -ex

# run_psfgen.sh

# generates topology files for pdb files in a directory.


if [ $# -ne 2 ]; then
    echo './run_psfgen.sh structure_file output_directory'
    exit 1
fi

module add vmd/1.9.3-fasrc01 readline

pdb=$1
output_directory=$2
: ${template:="/home01/home01/kibumpark/pkg/dbfold/MCPU/user_scripts/generate_topo.template.tcl"}

if [ ! -d "${output_directory}" ]; then
    mkdir -p "${output_directory}"
fi


echo "Working on ${pdb}"
fileroot=$(basename ${pdb} .pdb)
temptcl=$(mktemp -p ./)

sed \
    -e s:INPDBXXX:${pdb}: \
    -e s:OUTPSFXXX:"${output_directory}/${fileroot}.psf": \
    -e s:OUTPDBXXX:"${output_directory}/${fileroot}.pdb": \
    < "$template" \
    > $temptcl
    
cat $temptcl
vmd -dispdev text -e $temptcl
rm -f $temptcl
rm -f temp.pdb
