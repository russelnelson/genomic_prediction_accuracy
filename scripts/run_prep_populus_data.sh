#!/bin/bash -l
#SBATCH --mem=20G
#SBATCH -t 6-24:00:00
#SBATCH -J prep_populus

data=~/coalescence/prediction_accuracy/genomic_prediction_accuracy/data/populus_data/gatk_882_WG_genotypes_biallelic_snps_VQSR_non_filtered_bis.gt
phenos=~/coalescence/prediction_accuracy/genomic_prediction_accuracy/data/populus_data/gatk_882_WG_genotypes_biallelic_snps_VQSR_non_filtered_bis.gt.sorted_phenotypes.txt

sed '/0|1/ s//1/g' ${data}.vcf > ${data}.geno
sed -i '/0|0/ s//0/g' ${data}.geno
sed -i '/1|1/ s//2/g' ${data}.geno
sed -i '/1|0/ s//1/g' ${data}.geno

Rscript ~/coalescence/prediction_accuracy/genomic_prediction_accuracy/scripts/prep_populus_data.R ${data}.geno $phenos
