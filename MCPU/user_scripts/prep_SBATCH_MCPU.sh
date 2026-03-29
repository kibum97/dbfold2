#!/bin/bash
set -e

HELP="prep_SBATCH_MCPU.sh [options] input_protein output_directory
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

: ${MCPU_SHELL_TEMPLATE:="${MCPU_PATH}/user_scripts/shell_template"}

# Params for SBATCH
: ${USERVAR_NODES:=1}
: ${USERVAR_NTASKS:=16}
: ${USERVAR_JOB:="Default"}
: ${USERVAR_PARTITION:="shakhnovich,shared"} # and this should be the same
: ${USERVAR_MEM:=85000} #mem-per-cpu
: ${USERVAR_DAY:="1"}
: ${USERVAR_HOUR:="00"}
: ${USERVAR_MIN:="00"}
: ${USERVAR_DIR:="./"}

# input arguments
output_directory=$1

# some checks for files
if [ $# -lt 1 ]; then
    echo "not enough inputs, exiting"
    echo "required: output_directory"
    echo "$HELP"
    exit 10
fi

output_directory=$(readlink -f ${output_directory})
if [ ! -d "${output_directory}" ]; then
    echo mking dir ${output_directory}
    mkdir -p ${output_directory}
fi
# --------------------------------------------------

# Generate required files
# sh file
echo "Creating sh file"
sed \
-e "s:USERVAR_NODES:${USERVAR_NODES}:" \
-e "s:USERVAR_NTASKS:${USERVAR_NTASKS}:" \
-e "s:USERVAR_JOB:${USERVAR_JOB}:" \
-e "s:USERVAR_PARTITION:${USERVAR_PARTITION}:" \
-e "s:USERVAR_MEM:${USERVAR_MEM}:" \
-e "s:USERVAR_DAY:${USERVAR_DAY}:" \
-e "s:USERVAR_HOUR:${USERVAR_HOUR}:" \
-e "s:USERVAR_MIN:${USERVAR_MIN}:" \
-e "s:USERVAR_DIR:${USERVAR_DIR}:" \
< $MCPU_SHELL_TEMPLATE \
> "${output_directory}/MCPU.sh"
echo "Created protein-specific ${output_directory}/MCPU.sh"
