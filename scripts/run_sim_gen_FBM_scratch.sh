#!/bin/bash -l
#SBATCH -t 1:30:00
#SBATCH --output=sim_gen_%A_%a.out
#SBATCH -J sim_gen
#SBATCH --mem=20G
#SBATCH --array=1-10



geno_pheno=$1
outdir=$2
hmean=$3
hsd=$4
gpf=${geno_pheno##*/} # strip directory


module load gcc/9.2.0


# NOTE: MUST do this step first, so the next mkdir can be run without
# `-p'.  If mkdir is run with `-p' it always exits 0.
mkdir /scratch/$USER/ 2>/dev/null

dst="/scratch/$USER/$SLURM_ARRAY_JOB_ID"

HOST=$(hostname --short)

mkdir $dst
if [[ $? == 0 ]]; then
    # I am the first job to run
    echo "$SLURM_ARRAY_TASK_ID: making $dst, copying files"

    cp ${geno_pheno}.rds $dst
    cp ${geno_pheno}.subset.meta $dst
    cp ${geno_pheno}.subset.windows.RDS $dst
    cp ${geno_pheno}.subset.012.clean.transpose $dst
    cp ${geno_pheno}.subset.Gmat $dst
    cp ${geno_pheno}.subset.windows.RDS $dst
    cp ${geno_pheno}.rds $dst
    cp ${geno_pheno}.bk $dst
    cp ${geno_pheno}.geno.meta.RDS $dst


    echo $HOST >$dst/READY
else
    # I am not the first, wait for the first to copy all the files in
    while [[ ! -e $dst/READY ]]; do
        echo "$SLURM_ARRAY_TASK_ID: waiting for $dst/READY"
        # This sleep is required
        sleep 5
    done
    echo $HOST >>$dst/READY
fi

mkdir ${dst}/r_${SLURM_ARRAY_TASK_ID}
cd ${dst}/r_${SLURM_ARRAY_TASK_ID}

Rscript ~/coalescence/prediction_accuracy/genomic_prediction_accuracy/scripts/run_sim_gen_FBM.R ${dst}/${gpf}.rds ${outdir}/ABC_res.txt ${dst}/${gpf}.subset.meta  ${dst}/${gpf}.subset.windows.RDS ${dst}/${gpf}.subset.012.clean.transpose ${dst}/${gpf}.subset.Gmat $hmean $hsd ${outdir}/sim_gen_stats_r${SLURM_ARRAY_TASK_ID}.txt 


# Flag to indicate this array task is done.
echo $HOST >> $dst/DONE

# POST-all array job teardown here
# This sleep is important
sleep 60
if [[ $(grep -c $HOST $dst/DONE) == $(grep -c $HOST $dst/READY) ]]; then
    # I am the last array job to finish on this host.
    echo "$SLURM_ARRAY_TASK_ID: removing: $dst"
    rm -r $dst/
fi
