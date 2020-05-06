#!/bin/bash -l
#SBATCH --array=1-10000
#SBATCH --mem=12G
#SBATCH -t 2-24:00:00
#SBATCH -J esd_sim_test


infile=~/coalescence/prediction_accuracy/genomic_prediction_accuracy/ABC/ABC_input_fixed_number_flat_100_sites_scale_1.RDS
outname=fixed_number_flat_100_sites_scale_1_sim_gen
scheme_D_res=~/coalescence/prediction_accuracy/genomic_prediction_accuracy/ABC/fixed_flat_100sites_1scale/ABC_scheme_D_dual_opt_res_flat_scale1_sites100.RDS

module load gcc/9.2.0

mkdir r_${SLURM_ARRAY_TASK_ID}
cd r_${SLURM_ARRAY_TASK_ID}

Rscript ~/coalescence/prediction_accuracy/genomic_prediction_accuracy/ML/run_sim_gen.R $infile $scheme_D_res $outname $SLURM_ARRAY_TASK_ID
