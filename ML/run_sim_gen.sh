#!/bin/bash -l
#SBATCH -t 7-24:00:00
#SBATCH -J esd_sim_test
#SBATCH --mem=60G
#SBATCH -t 4-24:00:00
#SBATCH -n 24


infile=~/coalescence/prediction_accuracy/genomic_prediction_accuracy/ABC/ABC_input_scale_1_pi_9999_h_5_df_5.RDS
outname=ABC_input_scale_1_pi_9999_h_5_df_5.RDS
scheme_D_res=~/coalescence/prediction_accuracy/genomic_prediction_accuracy/ABC/pi_9999_scale_1_h_5_df_5/ABC_scheme_D_.RDS

module load gcc/9.2.0

Rscript ~/coalescence/prediction_accuracy/genomic_prediction_accuracy/ML/run_sim_gen.R $infile $scheme_D_res $outname
