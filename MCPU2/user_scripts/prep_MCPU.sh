#!/bin/bash
set -e

HELP="prep_MCPU.sh [options] input_protein output_directory
  Options: (much diminished here)
    -h : help
"
HELP2="
Input is a minimized protein structure with no hydrogen atoms. It
should exist as a pdb file with a single-chain structure as input. We
are almost ready to run MCPU.

There are several accompanying files that need to be generated to run
MCPU. This script generates those files.

1. FASTA file and .sec_str file
2. Potential-related files (by running save_triple)
3. cfg file

Script then copies relevant files to the run directory. Many
parameters for the cfg file can be set using environment variables:
i.e. execute this script with USERVAR_PARAM=VALUE preceding the
script (e.g. USERVAR_MC_STEPS=12345 ./prep_MCPU.sh).

@author: victor zhao yzhao01@g.harvard.edu
"

# fixed params:
MCPU_PATH=/n/home01/kibumpark/pkg/dbfold/MCPU/

SBATCH_SCRIPT="${MCPU_PATH}/user_scripts/MCPU.sh"
THREE_TO_ONE="${MCPU_PATH}/user_scripts/three_to_one.sed"

: ${MCPU_CFG_TEMPLATE:="${MCPU_PATH}/user_scripts/template.cfg"}

# Params for MCPU
: ${USERVAR_INDIR:="None"}
: ${USERVAR_MC_STEPS:=10000}
: ${USERVAR_MC_PRINT_STEPS:=2000}	    # this
: ${USERVAR_PRINT_PDB:=1}
: ${USERVAR_MC_PDB_PRINT_STEPS:=2000} # and this should be the same
: ${USERVAR_MC_REPLICA_STEPS:=20000}
: ${USERVAR_MAX_EXCHANGE:=20}
: ${USERVAR_MC_TEMP_MIN:=0.45}
: ${USERVAR_TEMP_STEP:=0.1}
: ${USERVAR_NODES_PER_TEMP:=1}
: ${USERVAR_NUMBER_OF_CONTACTS_MAX:=0}
: ${USERVAR_CONTACTS_STEP:=10}
: ${USERVAR_UMBRELLA:=1}
: ${USERVAR_BIAS_UMBRELLA:=0.02}
: ${USERVAR_CONSTRAINT_FILE:=""}
: ${USERVAR_CONSTRAINT_VALUE:=20}
: ${USERVAR_USE_CLUSTER:=0.0}
: ${USERVAR_MAX_CLUSTERSTEP:=0}
: ${USERVAR_CA_CUTOFF:=7}

NRUNS=1

# input arguments
while getopts ":N:h" opt; do
    case $opt in
        h) echo "$HELP" >&2; echo "${HELP2}" >&2; exit ;;
	N) echo "Number of runs: $OPTARG"; NRUNS=$OPTARG; ;;
        \?) echo "Invalid option: -$OPTARG" >&2; ;;
	:) echo "Invalid option: -$OPTARG requires an argument" 1>&2; exit 1 ;;
    esac
done
shift $((OPTIND-1))

input_protein=$1
output_directory=$2

# some checks for files
if [ $# -lt 2 ]; then
    echo "not enough inputs, exiting"
    echo "required: input_protein output_directory"
    echo "$HELP"
    exit 10
fi

if [ ! -f "${input_protein}" ]; then
    echo "Input protein file not found, exiting"
    echo "$HELP"
    exit 10
fi

input_protein="$(readlink -f "$input_protein")"
fileroot="$(basename "${input_protein}" .pdb)"

output_directory=$(readlink -f ${output_directory})
if [ ! -d "${output_directory}" ]; then
    echo mking dir ${output_directory}
    mkdir -p ${output_directory}
fi


# --------------------------------------------------

tempdir=$(mktemp -d)
trap "rm -rf $tempdir" 0 2 3 15	# some magic to remove tempdir at end
cd ${tempdir}

# Generate required files
# FASTA sequence from PDB
echo ">${fileroot}" > "${fileroot}.fasta"
grep ATOM < "${input_protein}" \
    | cut -c18-20,23-26 \
    | uniq \
    | awk '{print $1}' \
    | sed -f ${THREE_TO_ONE} \
    | tr -d '\n' \
    | fold \
	  >> "${fileroot}.fasta"
echo >> "${fileroot}.fasta"
echo "Created ${fileroot}.fasta"

# Generate .sec_str
n_residues=$(grep ATOM < "${input_protein}" \
		 | awk '{print $4}' \
		 | uniq \
		 | wc -l )
echo $(head -c $n_residues < /dev/zero | tr '\0' '0' ) \
     > "${fileroot}.sec_str"
echo $(head -c $n_residues < /dev/zero | tr '\0' 'C' ) \
     >> "${fileroot}.sec_str"
echo "Created ${fileroot}.sec_str"

# Run save_triple (slowest part of this script)
echo "Running save_triple"
ln -s ${MCPU_PATH}/mcpu_prep/{sct,triple}.energy ./
${MCPU_PATH}/mcpu_prep/save_triple ${fileroot}
echo "Created ${fileroot}.{triple,sctorsion}"

echo "moving files to output directory"
yes | mv ${fileroot}.{fasta,triple,sctorsion,sec_str} "${output_directory}"

cp -rv ${MCPU_PATH}/config_files "${output_directory}"
# Yes, we must copy this directory even though cfg points to data
#   files with the energy parameters b/c the program is hard-coded
#   to read a few files with path of "../config_files/FILE"


run_directory="${output_directory}/MCPU_run"
mkdir -p ${run_directory}	# snakemaek does create this dir I think
cp -v ${MCPU_PATH}/src_mpi_umbrella/fold_potential_mpi \
"${run_directory}"

# cfg file
echo "Creating cfg file"
sed \
-e "s:VAR_INPUT:${input_protein}:" \
-e "s:USERVAR_INDIR:${USERVAR_INDIR}:" \
-e "s:VAR_TEMPLATE:${output_directory}/nothing.template:" \
-e "s:VAR_OUTPUT:${run_directory}/${fileroot}:" \
-e "s:VAR_PARAM:${output_directory}/${fileroot}:" \
-e "s:USERVAR_MC_STEPS:${USERVAR_MC_STEPS}:" \
-e "s:USERVAR_MC_PRINT_STEPS:${USERVAR_MC_PRINT_STEPS}:" \
-e "s:USERVAR_PRINT_PDB:${USERVAR_PRINT_PDB}:" \
-e "s:USERVAR_MC_PDB_PRINT_STEPS:${USERVAR_MC_PDB_PRINT_STEPS}:" \
-e "s:USERVAR_MC_REPLICA_STEPS:${USERVAR_MC_REPLICA_STEPS}:" \
-e "s:USERVAR_MAX_EXCHANGE:${USERVAR_MAX_EXCHANGE}:" \
-e "s:USERVAR_MC_TEMP_MIN:${USERVAR_MC_TEMP_MIN}:" \
-e "s:USERVAR_TEMP_STEP:${USERVAR_TEMP_STEP}:" \
-e "s:USERVAR_NODES_PER_TEMP:${USERVAR_NODES_PER_TEMP}:" \
-e "s:USERVAR_UMBRELLA:${USERVAR_UMBRELLA}:" \
-e "s:USERVAR_BIAS_UMBRELLA:${USERVAR_BIAS_UMBRELLA}:" \
-e "s:USERVAR_NUMBER_OF_CONTACTS_MAX:${USERVAR_NUMBER_OF_CONTACTS_MAX}:" \
-e "s:USERVAR_CONTACTS_STEP:${USERVAR_CONTACTS_STEP}:" \
-e "s:USERVAR_CA_CUTOFF:${USERVAR_CA_CUTOFF}:" \
-e "s:USERVAR_CONSTRAINT_FILE:${USERVAR_CONSTRAINT_FILE}:" \
-e "s:USERVAR_CONSTRAINT_VALUE:${USERVAR_CONSTRAINT_VALUE}:" \
-e "s:USERVAR_USE_CLUSTER:${USERVAR_USE_CLUSTER}:" \
-e "s:USERVAR_MAX_CLUSTERSTEP:${USERVAR_MAX_CLUSTERSTEP}:" \
< $MCPU_CFG_TEMPLATE \
> "${run_directory}/cfg"
echo "Created protein-specific ${run_directory}/cfg"

rm -fv "${output_directory}/nothing.template"
touch "${output_directory}/nothing.template"
echo "Created ${output_directory}/nothing.template"
