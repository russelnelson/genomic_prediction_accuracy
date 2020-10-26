#!/bin/bash -l
#SBATCH -t 1:30:00
#SBATCH --mem=8G
#SBATCH -J ABC_FBM
#SBATCH --array=1-465


module load gcc/9.2.0

geno_pheno=$1
outdir=$2
hsd=$3
hmean=$4

Rscript ~/coalescence/prediction_accuracy/genomic_prediction_accuracy/scripts/run_ABC_FBM.R $geno_pheno ${outdir}/r${SLURM_ARRAY_TASK_ID}_ABC_out.txt $hsd $hmean

