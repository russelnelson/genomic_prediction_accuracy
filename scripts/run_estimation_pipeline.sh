#!/bin/bash -l

# ABC
j1=$(sbatch --parsable ~/coalescence/prediction_accuracy/genomic_prediction_accuracy/scripts/run_ABC.sh)

# gen sim
j2=$(sbatch  --dependency=afterok:$j1 --parsable ~/coalescence/prediction_accuracy/genomic_prediction_accuracy/scripts/run_gen_sim.sh)

# rf
j3=$(sbatch  --dependency=afterok:$j2 --parsable ~/coalescence/prediction_accuracy/genomic_prediction_accuracy/run_rf.sh)

# scale estimation
sbatch --parsable --dependency=afterok:$j3 ~/coalescence/prediction_accuracy/genomic_prediction_accuracy/run_regression.sh



