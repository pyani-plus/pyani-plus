#!/bin/bash

#=============================================================
#
# Job script for phylogenetic analysis
#
#=============================================================

#======================================================
# Propagate environment variables to the compute node
#SBATCH --export=ALL
# Project Account
#SBATCH --account=pritchard-grxiv
# Run in the standard partition (queue)
#SBATCH --partition=dev
# No of cores required (max. of 40)
#SBATCH --ntasks=12 --nodes=1
# Runtime (hard)
#SBATCH --time=00:30:00
# Job name
#SBATCH --job-name=snakemake_test
# Output file
#SBATCH --output=slurm-%j.out
#======================================================

# Enable conda enviroment
module load anaconda/python-3.8.8/2021.05
source ~/.bashrc

#Orthofinder only works in orthofinder conda environment
conda activate pyani_plus_smk
#=========================================================
# Prologue script to record job details
#=========================================================
/opt/software/scripts/job_prologue.sh
#----------------------------------------------------------

python -m pytest -v tests/snakemake/test_snakemake_slurm.py

conda deactivate

#=========================================================
# Epilogue script to record job endtime and runtime
#=========================================================
/opt/software/scripts/job_epilogue.sh
#----------------------------------------------------------
