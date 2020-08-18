#!/bin/bash -l
#SBATCH -t 1-12:00:00
#SBATCH -J ABC_FBM
#SBATCH --array=1-10000


module load gcc/9.2.0

geno_pheno=$1
outdir=$2
data_type=$3
hsd=$4
hmean=$5

Rscript ~/coalescence/prediction_accuracy/genomic_prediction_accuracy/scripts/run_ABC_FBM.R $geno_pheno ${outdir}/r${SLURM_ARRAY_TASK_ID}_ABC_out.txt $data_type $hsd $hmean

