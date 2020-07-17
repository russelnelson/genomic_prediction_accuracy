#!/bin/bash -l
#SBATCH --mem=128G
#SBATCH -t 4-24:00:00
#SBATCH -J ABC_imp_small_hsd
#SBATCH -n 24
#SBATCH -p bigmemm


module load gcc/9.2.0

geno_pheno=$1
outdir=$2
data_type=$3
hsd=$4

Rscript ~/coalescence/prediction_accuracy/genomic_prediction_accuracy/scripts/run_ABC.R $geno_pheno ${outdir}/ABC_out.RDS $data_type $hsd
