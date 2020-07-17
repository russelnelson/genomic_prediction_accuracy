#!/bin/bash -l

geno_pheno=~/coalescence/prediction_accuracy/genomic_prediction_accuracy/results/pi_.9999_scale_1_h_.5/ABC_input_scale_1_pi_9999_h_5_df_5.RDS
outdir=~/coalescence/prediction_accuracy/genomic_prediction_accuracy/results/pi_.9999_scale_1_h_.5/h_small_sd_imputed
data_type=imputed
hsd=0.02



# ABC
j1=$(sbatch --parsable ~/coalescence/prediction_accuracy/genomic_prediction_accuracy/scripts/run_ABC.sh $geno_pheno $outdir $data_type $hsd)

# gen sim
j2=$(sbatch --dependency=afterok:$j1 --parsable ~/coalescence/prediction_accuracy/genomic_prediction_accuracy/scripts/run_sim_gen.sh $geno_pheno $outdir $data_type)

# rf
j3=$(sbatch --dependency=afterok:$j2 --parsable ~/coalescence/prediction_accuracy/genomic_prediction_accuracy/scripts/run_rf.sh $geno_pheno $outdir $data_type)

# scale estimation
sbatch  --dependency=afterok:$j3 ~/coalescence/prediction_accuracy/genomic_prediction_accuracy/scripts/run_regression.sh $outdir



