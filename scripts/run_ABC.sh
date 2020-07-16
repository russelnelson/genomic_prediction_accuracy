#!/bin/bash -l
#SBATCH --array=1-1000
#SBATCH --mem=24G
#SBATCH -t 1-12:00:00
#SBATCH -J esd_ABC_test


infile=~/coalescence/prediction_accuracy/genomic_prediction_accuracy/data/theta4k_1000_10_rho40k.txt
outname=initial_test_ABC

module load gcc/9.2.0

mkdir r_${SLURM_ARRAY_TASK_ID}
cd r_${SLURM_ARRAY_TASK_ID}

Rscript ~/coalescence/prediction_accuracy/genomic_prediction_accuracy/run_ABC.R $infile $outname $SLURM_ARRAY_TASK_ID
