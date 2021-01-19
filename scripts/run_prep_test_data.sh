#!/bin/bash -l
#SBATCH --mem=4G
#SBATCH -t 4:00:00
#SBATCH -J prep_test
#SBATCH -p high

data=~/coalescence/prediction_accuracy/genomic_prediction_accuracy/data/test_data/bayesB_fixed_30_sites_h_.5_scale_1_trial_2.gt

sed '/0|1/ s//1/g' ${data}.vcf > ${data}.geno
sed -i '/0|0/ s//0/g' ${data}.geno
sed -i '/1|1/ s//2/g' ${data}.geno
sed -i '/1|0/ s//1/g' ${data}.geno

Rscript ~/coalescence/prediction_accuracy/genomic_prediction_accuracy/scripts/prep_test_data.R ${data}.geno  ${data}.phenos.txt
