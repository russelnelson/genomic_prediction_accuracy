#!/bin/bash -l
#SBATCH --array=1-10000
#SBATCH --mem=12G
#SBATCH -t 2-24:00:00
#SBATCH -J esd_sim_test


infile=~/coalescence/prediction_accuracy/genomic_prediction_accuracy/ABC/pi_.999_df_5_scale_1_h.75/ABC_input_.999pi_.75h_5df_scale1.RDS
outname=ABC_.999pi_.75h_5df_scale1_scale_not_fixed_sim_gen
scheme_D_res=~/coalescence/prediction_accuracy/genomic_prediction_accuracy/ABC/pi_.999_df_5_scale_1_h.75/ABC_scheme_D_dual_opt_res_pi.999_scale_1.RDS

module load gcc/9.2.0

mkdir r_${SLURM_ARRAY_TASK_ID}
cd r_${SLURM_ARRAY_TASK_ID}

Rscript ~/coalescence/prediction_accuracy/genomic_prediction_accuracy/ML/run_sim_gen.R $infile $scheme_D_res $outname $SLURM_ARRAY_TASK_ID
