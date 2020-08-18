#!/bin/bash -l
#SBATCH --mem=40G
#SBATCH -t 6-12:00:00
#SBATCH -n 24
#SBATCH -J make_GRM

gcta_exe=~/bin/gcta64 # path to the gcta executable
file_dir=~/coalescence/prediction_accuracy/genomic_prediction_accuracy/data/populus_data # the directory containing the plink files
pref=gatk_882_WG_genotypes_biallelic_snps_VQSR_non_filtered_bis.gt # plink file prefix
maf=0.05 # maf filter

${gcta_exe} --bfile ${file_dir}/${pref} --autosome-num 1528 --maf $maf --make-grm --out ${file_dir}/${pref} --thread-num 24
