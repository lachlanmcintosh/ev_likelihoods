#!/bin/bash
#SBATCH --job-name=singlecore_job_sage
#SBATCH --ntasks=1
#SBATCH --mem=2G
#SBATCH --requeue

# $1 refers to the first argument passed to the script
sage generate_all_path_combinations_for_given_length.sage "$1"

