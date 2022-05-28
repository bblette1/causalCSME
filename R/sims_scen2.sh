#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8g

srun R CMD BATCH --vanilla "--args seed$SLURM_ARRAY_TASK_ID" sims_scen2.R sims_scen2-$SLURM_ARRAY_TASK_ID.Rout
