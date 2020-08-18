#!/bin/bash -l
#SBATCH --array=1-1528
#SBATCH --mem=4G
#SBATCH -t 6-12:00:00
#SBATCH -J make_GRM

gcta_exe=~/bin/gcta64 # path to the gcta executable
list=~/coalescence/prediction_accuracy/genomic_prediction_accuracy/data/populus_data/chr_info.txt # A file containing the names of the chromosomes/scaffolds (need to be numbers) in column 2
file_dir=~/coalescence/prediction_accuracy/genomic_prediction_accuracy/data/populus_data # the directory containing the plink files
pref=gatk_882_WG_genotypes_biallelic_snps_VQSR_non_filtered_bis.gt # plink file prefix
maf=0.05 # maf filter


string="sed -n ${SLURM_ARRAY_TASK_ID}p ${list}" 
str=$($string)

var=$(echo $str | awk -F"\t" '{print $1, $2}')   
set -- $var
c1=$1
c2=$2

${gcta_exe}  --bfile ${file_dir}/${pref} --chr $c2 --autosome-num 1528 --maf $maf --make-grm --out ${file_dir}/chr${c2}_${pref}
