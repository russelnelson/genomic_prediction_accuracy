#!/bin/bash -l
#SBATCH --mem=40G
#SBATCH -t 4-24:00:00
#SBATCH -J pred
#SBATCH --array=1-100

Rscript run_gs.R /home/hemstrow/coalescence/prediction_accuracy/genomic_prediction_accuracy/results/test_data/sites_30_scale_1_h_.5/hsd.05/trial_1/gs_r${SLURM_ARRAY_TASK_ID}.txt
