#!/bin/bash -l

geno_pheno=~/coalescence/prediction_accuracy/genomic_prediction_accuracy/data/test_data/bayesB_fixed_30_sites_h_.5_scale_1_trial_2.gt
outdir=~/coalescence/prediction_accuracy/genomic_prediction_accuracy/results/test_data/sites_30_scale_1_h_.5/hsd.05/trial_2_non_subset
hsd=0.05
hmean=.5


# ABC
j1=$(sbatch --parsable ~/coalescence/prediction_accuracy/genomic_prediction_accuracy/scripts/run_ABC.sh $geno_pheno $outdir $data_type $hsd $hmean)

# gen sim
j2=$(sbatch --dependency=afterok:$j1 --parsable ~/coalescence/prediction_accuracy/genomic_prediction_accuracy/scripts/run_sim_gen.sh $geno_pheno $outdir $data_type)
#j2=$(sbatch --parsable ~/coalescence/prediction_accuracy/genomic_prediction_accuracy/scripts/run_sim_gen.sh $geno_pheno $outdir $data_type)


# rf
j3=$(sbatch --dependency=afterok:$j2 --parsable ~/coalescence/prediction_accuracy/genomic_prediction_accuracy/scripts/run_rf.sh $geno_pheno $outdir $data_type)
#j3=$(sbatch --parsable ~/coalescence/prediction_accuracy/genomic_prediction_accuracy/scripts/run_rf.sh $geno_pheno $outdir $data_type)

# scale estimation
sbatch  --dependency=afterok:$j3 ~/coalescence/prediction_accuracy/genomic_prediction_accuracy/scripts/run_regression.sh $outdir
#sbatch  ~/coalescence/prediction_accuracy/genomic_prediction_accuracy/scripts/run_regression.sh $outdir



