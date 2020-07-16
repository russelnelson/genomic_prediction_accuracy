#!/bin/bash -l
#SBATCH --mem=128G
#SBATCH -t 4-24:00:00
#SBATCH -J esd_ABC_test
#SBATCH -n 24
#SBATCH -p bigmemm


module load gcc/9.2.0

Rscript ~/coalescence/prediction_accuracy/genomic_prediction_accuracy/run_ABC_hyb_r1.R
