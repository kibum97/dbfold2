#!/bin/bash
#This script takes the last files of one simulation whose output path is specified in the first input variable, 
#and saves them as input files for another simulation in a directory specified in second input variable
#The MC timestep constituting the last files of the simulation is given by the third argument

#MCPU will always recognize files of the form 0.pdb, 1.pdb...as such, this reads myrank for each output file by
#looking at the log file, and then creates a new file myrank.pdb

#Sample instance would be . ./Starting_files_from_previous_simulation.sh 1igd/MultiUmbrella 1igd/MultiUmbrella 99500000

input_dir=$1
output_dir=$2
final_step=$3

array=($( ls $input_dir/*.log ))

for logfile in "${array[@]}"
do
	regex="^\.?/?(.+)_T_([0-9\.]+)_([0-9]+)\.log$"
    
    if [[ $logfile =~ $regex ]]; then
        protein="${BASH_REMATCH[1]}"
        temp="${BASH_REMATCH[2]}"
        setpoint="${BASH_REMATCH[3]}"
        
    else
        echo "Invalid filename format: $logfile"
        exit 1
    fi

	line=$(grep 'myrank' "$logfile")
    IFS=' ' read -r _ _ myrank <<< "$line"
    IFS=',' read -r myrank _ <<< "$myrank"

	cp $input_dir/${protein}_${temp}_$setpoint.$final_step $output_dir/$myrank.pdb
done