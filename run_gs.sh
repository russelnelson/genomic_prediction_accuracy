#!/bin/bash -l
#SBATCH --array=1-300
#SBATCH --mem=24G
#SBATCH -t 1-12:00:00
#SBATCH -J gsim


infile=theta4k_1000_10_rho40k.txt
outfile=outfiles/initial_runs.txt


Rscript prepare_for_growth.R $infile ${outfile}_${SLURM_ARRAY_TASK_ID} $SLURM_ARRAY_TASK_ID
