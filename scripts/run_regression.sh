#!/bin/bash -l
#SBATCH --mem=20G
#SBATCH -t 4-24:00:00
#SBATCH -J regression

module load gcc/9.2.0

outdir=$1

Rscript ~/coalescence/prediction_accuracy/genomic_prediction_accuracy/scripts/run_regression.R ${outdir}/RF.RDS ${outdir}/ABC_out.RDS ${outdir}/regression.RDS
