#!/bin/bash -l
#SBATCH --mem=128G
#SBATCH -t 4-24:00:00
#SBATCH -J prep_populus
#SBATCH -p bigmemm

Rscript ~/coalescence/prediction_accuracy/genomic_prediction_accuracy/scripts/prep_populus_data.R
