#!/bin/bash -l
#SBATCH -t 7-24:00:00
#SBATCH -J sim_gen
#SBATCH --mem=16G
#SBATCH -t 4-24:00:00
#SBATCH --array=1-10000



geno_pheno=$1
outdir=$2
data_type=$3

module load gcc/9.2.0

mkdir ${outdir}/r_${SLURM_ARRAY_TASK_ID}
cd ${outdir}/r_${SLURM_ARRAY_TASK_ID}

Rscript ~/coalescence/prediction_accuracy/genomic_prediction_accuracy/scripts/run_sim_gen.R $geno_pheno ${outdir}/ABC_out.RDS ${outdir}/sim_gen_stats_r${SLURM_ARRAY_TASK_ID}.txt $data_type

cd ~/
rm -r ${outdir}/r_${SLURM_ARRAY_TASK_ID}
