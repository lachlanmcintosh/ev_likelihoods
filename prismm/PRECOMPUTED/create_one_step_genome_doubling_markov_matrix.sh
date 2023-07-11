#!/bin/bash
#SBATCH --job-name=create_G
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00  # Adjust as needed

sage create_one_step_genome_doubling_markov_matrix.sage "$1"  # Assuming max_CN is passed as an argument

