#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 48
#SBATCH -p test
#SBATCH -t 0-12:00
#SBATCH --mem=180G
#SBATCH -J landscape                                                            
#SBATCH -o logs/landscape_%j.out                                                         
#SBATCH -e logs/landscape_%j.err

set -e

module load Mambaforge

mamba activate dbfold2

cd /n/home01/kibumpark/group_folder/p19.7__Actin_linker_simulations/ACTA

python /n/home01/kibumpark/group_folder/p19.7__Actin_linker_simulations/scripts/compute_mbar.py \
--pdbroots acta_l \
--eq_step 600000000 \
--basedir /n/home01/kibumpark/group_folder/p19.7__Actin_linker_simulations/ACTA \
--homedir /n/home01/kibumpark/group_folder/p19.7__Actin_linker_simulations/ACTA/Analysis/ \
--cache

python /n/home01/kibumpark/group_folder/p19.7__Actin_linker_simulations/scripts/compute_mbar.py \
--pdbroots actb_l \
--eq_step 600000000 \
--basedir /n/home01/kibumpark/group_folder/p19.7__Actin_linker_simulations/ACTB \
--homedir /n/home01/kibumpark/group_folder/p19.7__Actin_linker_simulations/ACTB/Analysis/ \
--cache