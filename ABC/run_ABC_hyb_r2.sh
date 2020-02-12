#!/bin/bash -l
#SBATCH --array=1-2500
#SBATCH --mem=12G
#SBATCH -t 2-24:00:00
#SBATCH -J esd_ABC_test


infile=~/coalescence/prediction_accuracy/genomic_prediction_accuracy/ABC/pi_.9999_df_5_scale_1_h.75/ABC_test_pi.9999_df_5_scale_1_h_.75.RDS
outname=ABC_.9999pi_.75h_5df_scale1_scale_not_fixed_hyb_out
scheme_D_res=~/coalescence/prediction_accuracy/genomic_prediction_accuracy/ABC/pi_.9999_df_5_scale_1_h.75/ABC_scheme_D_dual_opt_res_pi.9999_scale_1.RDS

module load gcc/9.2.0

mkdir r_${SLURM_ARRAY_TASK_ID}
cd r_${SLURM_ARRAY_TASK_ID}

Rscript ~/coalescence/prediction_accuracy/genomic_prediction_accuracy/run_ABC_hyb_r2.R $infile $scheme_D_res $outname $SLURM_ARRAY_TASK_ID
