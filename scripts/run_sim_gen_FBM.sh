#!/bin/bash -l
#SBATCH -t 7-24:00:00
#SBATCH -J sim_gen
#SBATCH --mem=2G
#SBATCH -t 4-24:00:00
#SBATCH --array=1-10000



geno_pheno=$1
outdir=$2
hmean=$3
hsd=$4

module load gcc/9.2.0

mkdir ${outdir}/r_${SLURM_ARRAY_TASK_ID}
cd ${outdir}/r_${SLURM_ARRAY_TASK_ID}

Rscript ~/coalescence/prediction_accuracy/genomic_prediction_accuracy/scripts/run_sim_gen_FBM.R ${geno_pheno}.rds ${outdir}/ABC_res.txt ${geno_pheno}.subset.meta  ${geno_pheno}.subset.windows.RDS ${geno_pheno}.subset.012.clean.transpose ${geno_pheno}.subset.Gmat $hmean $hsd ${outdir}/sim_gen_stats_r${SLURM_ARRAY_TASK_ID}.txt 


cd ~/
rm -r ${outdir}/r_${SLURM_ARRAY_TASK_ID}
