#!/bin/bash

# Define the maximum length of the path
max_path_length=8

# Define the maximum number of jobs
max_jobs=400

# Generate all paths that we want to compute
srun generate_all_path_combinations_for_given_length.sh "$max_path_length"

# Generate the one step anueploidy matrix that will be used to compute path likelihoods
srun create_one_step_anueploidy_markov_matrix.sh "$max_path_length"

# Generate the genome doubling matrix
srun create_one_step_genome_doubling_markov_matrix.sh "$max_path_length"

# Calculate the maximum copy number
max_CN=$((2 ** $max_path_length))
path_description="p$max_path_length_v4"

submit_job() {
  local p_up=$1
  local p_down=$2

  FILE="MATRICES/subbed_u$p_up_d$p_down_$path_description.csv"

  if ! [ -f "$FILE" ]; then
    echo "$FILE doesn't exist!"

    local num_jobs
    num_jobs=$(squeue -u "$USER" | awk 'NR > 1' | wc -l)

    while [ $num_jobs -ge $max_jobs ]; do
      sleep 10
      num_jobs=$(squeue -u "$USER" | awk 'NR > 1' | wc -l)
    done

    sbatch ./generate_path_likelihoods.sh "$p_up" "$p_down" "$max_CN" "$max_path_length" "$path_description"
  fi
}

for rep in {0..0}; do
  for p_up in {0..100}; do
    for p_down in {0..100}; do
      if [ $((p_up+p_down)) -lt 101 ]; then
        submit_job "$p_up" "$p_down"
      fi
    done
  done
done

# Check if all files exist
python check_all_files_exist.py

