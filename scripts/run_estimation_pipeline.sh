#!/bin/bash -l

# prep for non-FBM data: run make_test_data.R with correct settings, unzip the resulting vcf, run prep_test_data.R, then run convert_FBM_to_rds.R
# This will create all of the needed datasets!
# Lastly, make sure that the outdir directory listed here actually exists

geno_pheno=/home/hemstrow/coalescence/prediction_accuracy/genomic_prediction_accuracy/data/test_data/bayesB_fixed_30_sites_h_.5_scale_1_trial_3_genopheno.RDS
outdir=~/coalescence/prediction_accuracy/genomic_prediction_accuracy/results/test_data/sites_30_scale_1_h_.5/hsd.05/t1_re/
hsd=0.05
hmean=.5
peak_delta=0.05
peak_pcut=0.001

# ABC
j1=$(sbatch --parsable ~/coalescence/prediction_accuracy/genomic_prediction_accuracy/scripts/run_ABC.sh $geno_pheno $outdir $hsd $hmean)

# gen sim
j2=$(sbatch --dependency=afterok:$j1 --parsable ~/coalescence/prediction_accuracy/genomic_prediction_accuracy/scripts/run_sim_gen.sh $geno_pheno $outdir $hsd $hmean $peak_delta $peak_pcut)
#j2=$(sbatch --parsable ~/coalescence/prediction_accuracy/genomic_prediction_accuracy/scripts/run_sim_gen.sh $geno_pheno $outdir $hsd $hmean $peak_delta $peak_pcut)

# rf
j3=$(sbatch --dependency=afterok:$j2 --parsable ~/coalescence/prediction_accuracy/genomic_prediction_accuracy/scripts/run_rf.sh $geno_pheno $outdir $peak_delta $peak_pcut)
#j3=$(sbatch --parsable ~/coalescence/prediction_accuracy/genomic_prediction_accuracy/scripts/run_rf.sh $geno_pheno $outdir $peak_delta $peak_pcut)

# scale estimation
sbatch  --dependency=afterok:$j3 ~/coalescence/prediction_accuracy/genomic_prediction_accuracy/scripts/run_regression.sh $outdir
#sbatch  ~/coalescence/prediction_accuracy/genomic_prediction_accuracy/scripts/run_regression.sh $outdir



