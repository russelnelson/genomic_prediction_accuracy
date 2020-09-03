#!/bin/bash -l
#SBATCH -t 1:00:00
#SBATCH -J cat_abc
#SBATCH -p high

outdir=$1
# cat together the ABC res!

cat ${outdir}/r*_ABC_out.txt > ${outdir}/ABC_res.txt


