#!/bin/bash -l
#SBATCH -t 7-24:00:00
#SBATCH -J sim_gen
#SBATCH --mem=10G
#SBATCH -t 4-24:00:00
#SBATCH --array=1-100000



geno_pheno=$1
outdir=$2
data_type=$3

module load gcc/9.2.0

mkdir r_${SLURM_ARRAY_TASK_ID}
mv r_${SLURM_ARRAY_TASK_ID}

Rscript ~/coalescence/prediction_accuracy/genomic_prediction_accuracy/ML/run_sim_gen.R $geno_pheno ${outdir}/ABC_out.RDS ${outdir}/sim_gen_stats_r${SLURM_ARRAY_TASK_ID}.txt $data_type

mv ..
rm -r r_${SLURM_ARRAY_TASK_ID}