#!/bin/bash
set -ex

<<'EOF'

After NAMD 5000 step minimization on the structures,
1. Create from template a TCL script to extract minimized structure
2. Call vmd to run script
3. Take output temp pdb and run it through sed to change some aliased names

EOF

echo extract_minimized_structure.sh; echo

if [ $# -ne 3 ]; then
    echo "./extract_minimized_structure.sh structure trajectory output_dir"
    exit 1
fi

module purge
module add vmd/1.9.3-fasrc01 readline

structure=$1
trajectory=$2
output_directory=$3
template=/n/home00/vzhao/code/dbfold/MCPU/user_scripts/extract_pdb.template.tcl

if [ ! -d "${output_directory}" ]; then
    mkdir -p "${output_directory}"
fi

# rm -rf MIN_OUT_PDB_TEMP.pdb extract_minimized_structure.tcl
fileroot=$(basename ${structure} .psf)
echo "Working on ${fileroot}"

outfile=${output_directory}/${fileroot}.pdb

temptclscript=$(mktemp -p ./)
echo "temporary tcl to $temptclscript"

sed \
    -e s:INPSFXXX:"${structure}": \
    -e s:INDCDXXX:"${trajectory}": \
    -e s:OUTPDBXXX:"${outfile}": \
    < "$template" \
    > $temptclscript
if [ ! $? ]; then
    echo "sed failed?"
    exit 1
fi

vmd -dispdev text -e $temptclscript

if [ $? ]; then
    sed -i \
	-e s/HSE/HIS/ \
	-e "/ILE/ s/CD /CD1/" \
	${outfile}
else
    echo "extract must have failed"
    exit 1
fi

rm -f $temptclscript
