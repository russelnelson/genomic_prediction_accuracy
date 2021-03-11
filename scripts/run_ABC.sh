#!/bin/bash -l
#SBATCH --mem=60G
#SBATCH -t 6-24:00:00
#SBATCH -J ABC
#SBATCH -n 24
#SBATCH -p bigmemh

module load gcc/9.2.0

geno_pheno=$1
outdir=$2
hsd=$3
hmean=$4

Rscript ~/coalescence/prediction_accuracy/genomic_prediction_accuracy/scripts/run_ABC.R $geno_pheno ${outdir}/ABC_res.txt $hsd $hmean

