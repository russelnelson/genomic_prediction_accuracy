#!/bin/bash -l
#SBATCH --mem=128G
#SBATCH -t 4-24:00:00
#SBATCH -J rf
#SBATCH -n 24
#SBATCH -p bigmemm


module load gcc/9.2.0

geno_pheno=$1
outdir=$2
data_type=$3

cat ${outdir}/sim_gen_stats_r*.txt > ${outdir}/sim_gen_stats.txt
rm ${outdir}/sim_gen_stats_r*.txt

Rscript ~/coalescence/prediction_accuracy/genomic_prediction_accuracy/scripts/run_rf.R ${outdir}/sim_gen_stats.txt $geno_pheno ${outdir}/RF.RDS $data_type
