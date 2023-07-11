#!/bin/bash

#SBATCH --job-name=simulation_job
#SBATCH --output=simulation_job.out
#SBATCH --error=simulation_job.err
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=100GB
#SBATCH --nodes=1
#SBATCH --partition=regular

# Change these values as per your need
BEST_EST_CUTOFF1=3.1
BEST_EST_CUTOFF2=7.8
ARB_CUTOFF1=3.8
ARB_CUTOFF2=7.6


srun python simulation_sets.py --base_filename simgrid \
    --cutoffs_best_estimate ${BEST_EST_CUTOFF1} ${BEST_EST_CUTOFF2} \
    --cutoffs_arbitrary_estimate ${ARB_CUTOFF1} ${ARB_CUTOFF2}

