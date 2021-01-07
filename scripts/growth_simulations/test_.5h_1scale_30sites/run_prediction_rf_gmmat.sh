#!/bin/bash -l
#SBATCH --mem=128G
#SBATCH -t 4-24:00:00
#SBATCH -J pred
#SBATCH -n 24
#SBATCH -p bigmemm


Rscript run_prediction_rf_gmmat.R
