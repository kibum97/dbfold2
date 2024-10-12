#!/bin/bash
set -ex

HELP="run_MCPU-multiconf.sh [options] input_dir output_directory
  Options:
    -h : help
    -n NPROC : specify number of processors and indirectly
               specify the size of the temperature ladder
    -N NRUNS : specify number of runs to do
    -l LEN   : number of steps in a run
    -e EVERY : save structures for every specified trajectories
               give non-number as argument for no saving at all.
    -o FREQ  : save structures every FREQ steps
    -m       : save min E and min RMSD PDBs at end of sim
    -s       : skip preparation steps and do job submission
"
# ./run_MCPU.sh structure output_directory
echo run_MCPU-multiconf.sh; echo
<<"EOF"
After preparing minimized protein structures with hydrogen atoms
removed, we are ready to run MCPU.

Based off of the readme.txt instructions as well as
produce_para_remc.pl

1. Requires an input protein as a .pdb file. no hydrogens.
   Single chain!
2. Accompanying input files will be generated in the directory
   containing the input protein.
3. Output will be placed in output_directory, with files named
   prefix.*, where prefix is the input pdb file minus .pdb suffix
# 4. Only one instance of this script should be running. Prep script
#    save_triple runs in the source dir, and file cfg in the source dir
#    also gets modified. Wouldn't want to modify those files while
#    things are running.

Item 4 has been corrected.

The number of processors dictates the temperature range of the
simulation because delta T is 0.1, with the
starting temperature being 0.1. So -n 32 is 0.1 to 3.2

@author: victor zhao yzhao01@g.harvard.edu
EOF

# fixed params:
MCPU_PATH=/n/home00/vzhao/code/dbfold/MCPU

THREE_TO_ONE="${MCPU_PATH}/user_scripts/three_to_one.sed"
SBATCH_SCRIPT="MCPU.sh"
MCPU_CFG_TEMPLATE="${MCPU_PATH}/src_mpi/template.cfg"
# within MCPU, files src_mpi/TEMPLATE.cfg has been edited
#  to accomodate this script.

: ${PARTITION:=shakhnovich,shakgpu,shared}
: ${MEMPERCPU:=1000}
: ${ALLOCTIME:=240}
: ${NPROC:=32}
: ${LENGTH:=2000000}
: ${EVERY:=10}                        # Every N runs will save structure trajectories
: ${OUTFREQ:=100000}                  # print pdb file freq
: ${DRYRUN:=}
SAVEMIN=0

# input arguments
while getopts ":hmsn:N:l:e:o:" opt; do
    case $opt in
        h) echo "$HELP" >&2; exit ;;
	m) echo "Save min E and min RMSD pdbs"; SAVEMIN=1 ;;
	s) echo "Skip prep"; SKIPPREP=TRUE ;;
        e) echo "Save every $OPTARG trajectories"; EVERY=$OPTARG; ;;
        n) echo "Using $OPTARG processors"; NPROC=$OPTARG; ;;
        # N) echo "Number of runs: $OPTARG"; NRUNS=$OPTARG; ;;
        l) echo "Length of run: $OPTARG"; LENGTH=$OPTARG; ;;
        o) echo "Output pdb every $OPTARG steps"; OUTFREQ=$OPTARG; ;;
        \?) echo "Invalid option: -$OPTARG" >&2; ;;
    esac
done
shift $((OPTIND-1))

input_directory=$1
output_directory=$2

# some checks for files
if [ $# -lt 2 ]; then
    echo "not enough inputs, exiting"
    echo "required: input_protein output_directory"
    echo "$HELP"
    exit 10
fi

if [ ! -d "${input_directory}" ]; then
    echo "Input directory does not exist, exiting"
    echo "$HELP"
    exit 10
fi

input_directory="$(readlink -f "$input_directory")"
file0="$(ls ${input_directory}/*.pdb | head -n1)"

if [ ! -f "$file0" ]; then
    echo "No pdbfiles in input dir, exiting"
    echo "$HELP"
    exit 10
fi

fileroot="$(basename "${input_directory}")"
# input_prefix="${input_directory}/${fileroot}"

if [ ! -d "${output_directory}" ]; then
    mkdir -p ${output_directory}
fi
output_directory=$(readlink -f ${output_directory})


# Begin --------------------------------------------------

if [ -z "${SKIPPREP}" ]; then
    tempdir=$(mktemp -d)
    trap "rm -rf $tempdir" 0 2 3 15
    echo "change directory: ${tempdir}"; cd ${tempdir}

    # Generate required files
    # FASTA sequence from PDB
    echo ">${fileroot}" > "${fileroot}.fasta"
    grep ATOM < "${file0}" \
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
    n_residues=$(grep ATOM < "${file0}" \
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
    # cp "${input_prefix}.fasta" ${tempdir}
    # echo "change directory: ${tempdir}"; cd ${tempdir}
    ln -s ${MCPU_PATH}/mcpu_prep/{sct,triple}.energy ./
    # ln -s ${MCPU_PATH}/mcpu_prep/save_triple ./
    ${MCPU_PATH}/mcpu_prep/save_triple ${fileroot}
    yes | mv ${fileroot}.{fasta,triple,sctorsion,sec_str} "${input_directory}"
    echo "Created ${fileroot}.{triple,sctorsion}"

    # echo "moving files to output directory"
    # cp -rv ${MCPU_PATH}/config_files "${output_directory}"
    # cp -v ${MCPU_PATH}/src_mpi/fold_potential_mpi "${output_directory}"
    # cp -v "${MCPU_PATH}/user_scripts/${SBATCH_SCRIPT}" "${output_directory}"
fi

echo "moving files to output directory"
cp -rv ${MCPU_PATH}/config_files "${output_directory}"
cp -v ${MCPU_PATH}/src_mpi/fold_potential_mpi "${output_directory}"
cp -v "${MCPU_PATH}/user_scripts/${SBATCH_SCRIPT}" "${output_directory}"

i=0

# For each run
for infile in $(ls -v ${input_directory}/*.pdb); do
    input_file="$(readlink -f "$infile")"

    # make subdirectory
    subdirectory="${output_directory}/run_${i}"
    mkdir "${subdirectory}"
    
    PRINT_PDB=0
    # this part in a subshell in case EVERY is non-number
    if ! [ "$EVERY" -eq "$EVERY" ] 2> /dev/null; then
	PRINT_PDB=0
    elif [ $(($i % $EVERY)) -eq 0 ]; then
	PRINT_PDB=1
    else
        PRINT_PDB=0
    fi

    # configure files
    sed \
	-e "s:VAR_OUTPUT:${subdirectory}/${fileroot}:" \
	-e "s:VAR_TEMPLATE:${subdirectory}:" \
	-e "s:VAR_INPUT:${input_file}:" \
	-e "s:VAR_PARAM:${input_directory}/${fileroot}:" \
	-e "s:VAR_STEPS:${LENGTH}:" \
	-e "s:VAR_PRINT:${PRINT_PDB}:" \
	-e "s:VAR_OUTFREQ:${OUTFREQ}:" \
	-e "s:VAR_SAVEMIN:${SAVEMIN}:" \
	< $MCPU_CFG_TEMPLATE \
	> "${subdirectory}/cfg"
    echo "Created protein-specific ${subdirectory}/cfg"

    rm -fv "${subdirectory}/nothing.template"
    touch "${subdirectory}/nothing.template"
    echo "Created ${subdirectory}/nothing.template"

    echo "change directory: $subdirectory"; cd $subdirectory
	
    echo "Now run program."
    # $([ -n "$DRYRUN" ] && echo echo) \
    # 	sbatch -p $PARTITION -n $NPROC -t $ALLOCTIME --mem-per-cpu $MEMPERCPU \
    # 	-J "MCPU_${fileroot}_R${i}" \
    # 	-o "MCPU_${fileroot}_run${i}.slurm" \
    # 	--wrap "mpirun -n $NPROC ./fold_potential_mpi cfg"
    $([ -n "$DRYRUN" ] && echo echo) \
	sbatch -p $PARTITION -n $NPROC -t $ALLOCTIME --mem-per-cpu $MEMPERCPU \
    	-J "MCPU_${fileroot}_R${i}" -o "MCPU_${fileroot}_run${i}.slurm" \
	../${SBATCH_SCRIPT} cfg

    echo "Job submitted, output: ${subdirectory}"
    let ++i
    sleep 0.5
done

<<NOTES

What's needed is an "installation" of MCPU in which everything is
configured properly. 

0. Put it in a fixed path location.
1. Change /PATHNAME/ instances to the path in config files
2. Set output directory
3. Edit configuration files to taste
4. I think some configuration settings have to be left at run time?
5. define.h parameters probably are fine?

NOTES
