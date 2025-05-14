#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 48
#SBATCH -p sapphire,test
#SBATCH -t 0-12:00
#SBATCH --mem=500G
#SBATCH -J landscape                                                            
#SBATCH -o logs/landscape_%j.out                                                         
#SBATCH -e logs/landscape_%j.err

set -e

module load Mambaforge

mamba activate dbfold2

cd /n/home01/kibumpark/group_folder/p19__Actin/Analysis/

python /n/home01/kibumpark/group_folder/p19__Actin/dbfold2/dbfold/scripts/compute_dbfold_landscape.py \
--pdbroots ACTB_cryo \
--eq_step 600000000 \
--substrcutre_dict /n/home01/kibumpark/group_folder/p19__Actin/Analysis/substructures/ACTB_substructures_definition_v2.pkl \
--basedir /n/home01/kibumpark/group_folder/p19__Actin/ \
--homedir /n/home01/kibumpark/group_folder/p19__Actin/Analysis/ \
--savedir /n/home01/kibumpark/group_folder/p19__Actin/Analysis/finer_def_figures/ \
--cache

python /n/home01/kibumpark/group_folder/p19__Actin/dbfold2/dbfold/scripts/compute_dbfold_landscape.py \
--pdbroots ACTA_cryo \
--eq_step 600000000 \
--substrcutre_dict /n/home01/kibumpark/group_folder/p19__Actin/Analysis/substructures/ACTA_substructures_definition_v2.pkl \
--basedir /n/home01/kibumpark/group_folder/p19__Actin/ \
--homedir /n/home01/kibumpark/group_folder/p19__Actin/Analysis/ \
--savedir /n/home01/kibumpark/group_folder/p19__Actin/Analysis/finer_def_figures/
--cache