#!/bin/bash -l
#SBATCH --mem=200G
#SBATCH -t 2-24:00:00
#SBATCH -J rf
#SBATCH -n 24
#SBATCH -p bigmemm


module load gcc/9.2.0

geno_pheno=$1
outdir=$2

cat ${outdir}/sim_gen_stats_r*.txt >> ${outdir}/sim_gen_stats.txt
#rm ${outdir}/sim_gen_stats_r*.txt

Rscript ~/coalescence/prediction_accuracy/genomic_prediction_accuracy/scripts/run_rf_FBM.R ${outdir}/sim_gen_stats.txt ${geno_pheno}.subset.meta ${geno_pheno}.subset.windows.RDS ${geno_pheno}.subset.012.clean.transpose ${geno_pheno}.subset.Gmat ${geno_pheno}.geno.meta.RDS ${outdir}/RF.RDS 
