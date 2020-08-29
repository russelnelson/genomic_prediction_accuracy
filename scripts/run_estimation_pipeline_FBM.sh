#!/bin/bash -l

geno_pheno=~/coalescence/prediction_accuracy/genomic_prediction_accuracy/data/populus_data/gatk_882_WG_genotypes_biallelic_snps_VQSR_non_filtered_bis.gt
outdir=~/coalescence/prediction_accuracy/genomic_prediction_accuracy/results/populus/h_0.05_sd
data_type=imputed
hsd=0.05
hmean=.5
sub_maf=.05
sub_thin=1000

# ABC
#j1=$(sbatch --parsable ~/coalescence/prediction_accuracy/genomic_prediction_accuracy/scripts/run_ABC_FBM.sh $geno_pheno $outdir $data_type $hsd $hmean)

# subset
#j2=$(sbatch --parsable prep_populus_GWAS_subset.sh $geno_pheno $sub_thin $sub_maf)

# gen sim
#j3=$(sbatch --dependency=afterok:$j2 --parsable ~/coalescence/prediction_accuracy/genomic_prediction_accuracy/scripts/run_sim_gen.sh $geno_pheno $outdir $data_type)
j3=$(sbatch --parsable ~/coalescence/prediction_accuracy/genomic_prediction_accuracy/scripts/run_sim_gen_FBM.sh $geno_pheno $outdir $hmean $hsd)


# rf
#j4=$(sbatch --dependency=afterok:$j3 --parsable ~/coalescence/prediction_accuracy/genomic_prediction_accuracy/scripts/run_rf.sh $geno_pheno $outdir $data_type)
#j4=$(sbatch --parsable ~/coalescence/prediction_accuracy/genomic_prediction_accuracy/scripts/run_rf.sh $geno_pheno $outdir $data_type)

# scale estimation
#sbatch  --dependency=afterok:$j4 ~/coalescence/prediction_accuracy/genomic_prediction_accuracy/scripts/run_regression.sh $outdir
#sbatch  ~/coalescence/prediction_accuracy/genomic_prediction_accuracy/scripts/run_regression.sh $outdir



