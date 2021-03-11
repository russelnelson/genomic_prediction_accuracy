#!/bin/bash -l
#SBATCH -t 7-24:00:00
#SBATCH -J sim_gen
#SBATCH --mem=24G
#SBATCH -t 4-24:00:00
#SBATCH --array=1-20000



geno_pheno=$1
outdir=$2
hsd=$3
hmean=$4
peak_delta=$5
peak_pcut=$6

module load gcc/9.2.0

mkdir ${outdir}/r_${SLURM_ARRAY_TASK_ID}
cd ${outdir}/r_${SLURM_ARRAY_TASK_ID}

Rscript ~/coalescence/prediction_accuracy/genomic_prediction_accuracy/scripts/run_sim_gen.R $geno_pheno ${outdir}/ABC_res.txt ${outdir}/sim_gen_stats_r${SLURM_ARRAY_TASK_ID}.txt $hsd $hmean $peak_delta $peak_pcut

cd ~/
rm -r ${outdir}/r_${SLURM_ARRAY_TASK_ID}
