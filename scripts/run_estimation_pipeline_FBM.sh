#!/bin/bash -l

geno_pheno=~/coalescence/prediction_accuracy/genomic_prediction_accuracy/data/test_data/bayesB_fixed_30_sites_h_.5_scale_1_trial_2.gt
outdir=~/coalescence/prediction_accuracy/genomic_prediction_accuracy/results/test_data/sites_30_scale_1_h_.5/hsd.05/trial_2
hsd=0.05
hmean=.5
sub_maf=.05
sub_thin=500
datamash_path=~/bin/datamash

# ABC and cat the results
#j1=$(sbatch --parsable ~/coalescence/prediction_accuracy/genomic_prediction_accuracy/scripts/run_ABC_FBM_scratch.sh $geno_pheno $outdir $hsd $hmean)
#j1c=$(sbatch --parsable  --dependency=afterok:$j1 ~/coalescence/prediction_accuracy/genomic_prediction_accuracy/scripts/cat_ABC.sh $outdir)

# subset
#j2=$(sbatch --parsable prep_GWAS_subset.sh $geno_pheno $sub_thin $sub_maf $outdir $datamash_path)

# gen sim
#j3=$(sbatch --parsable --dependency=afterok:$j2:$j1c ~/coalescence/prediction_accuracy/genomic_prediction_accuracy/scripts/run_sim_gen_FBM_scratch.sh $geno_pheno $outdir $hmean $hsd)
#j3=$(sbatch --parsable ~/coalescence/prediction_accuracy/genomic_prediction_accuracy/scripts/run_sim_gen_FBM_scratch.sh $geno_pheno $outdir $hmean $hsd)

# rf
j4=$(sbatch --parsable ~/coalescence/prediction_accuracy/genomic_prediction_accuracy/scripts/run_rf_FBM.sh $geno_pheno $outdir)
#j4=$(sbatch --dependency=afterok:$j3 --parsable ~/coalescence/prediction_accuracy/genomic_prediction_accuracy/scripts/run_rf_FBM.sh $geno_pheno $outdir)

# scale estimation
sbatch  --dependency=afterok:$j4 ~/coalescence/prediction_accuracy/genomic_prediction_accuracy/scripts/run_regression.sh $outdir
#sbatch ~/coalescence/prediction_accuracy/genomic_prediction_accuracy/scripts/run_regression.sh $outdir



