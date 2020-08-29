#!/bin/bash -l
#SBATCH --mem=40G
#SBATCH -t 4-24:00:00
#SBATCH -J prep_GWAS_populus
#SBATCH -p bigmemm

module load vcftools

data=~/coalescence/prediction_accuracy/genomic_prediction_accuracy/data/populus_data/gatk_882_WG_genotypes_biallelic_snps_VQSR_non_filtered_bis.gt
data=$1
thin=$2
maf=$3

#vcftools --gzvcf ${data}.vcf.gz --out ${data}.subset --maf 0.05 --thin 1000 --recode --recode-INFO-all

#vcftools --vcf ${data}.subset.recode.vcf --out ${data}.subset --012

#cut -f 2- ${data}.subset.012 > ${data}.subset.012.clean

#Rscript make_G_agh_matrix.R ${data}.subset.012.clean ${data}.subset.Gmat

#cut -f 1-2 ${data}.subset.recode.vcf | grep -v "#" - > ${data}.subset.meta

Rscript save_populus_windows.R ${data}.subset.meta ${data}.subset.windows.RDS


