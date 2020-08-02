#!/bin/bash -l
#SBATCH --mem=200G
#SBATCH -t 4-24:00:00
#SBATCH -J prep_populus
#SBATCH -p bigmemm

data=~/coalescence/prediction_accuracy/genomic_prediction_accuracy/data/populus_data/gatk_882_WG_genotypes_biallelic_snps_VQSR_non_filtered_bis.gt.vcf

sed -i '/0\\|1/ s//1/g' $data
sed -i '/0\\|0/ s//0/g' $data
sed -i '/1\\|1/ s//2/g' $data
sed -i '/1\\|0/ s//1/g' $data

Rscript ~/coalescence/prediction_accuracy/genomic_prediction_accuracy/scripts/prep_populus_data.R $data
