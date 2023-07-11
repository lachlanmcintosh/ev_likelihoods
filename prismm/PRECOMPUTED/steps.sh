#!/bin/bash

max_path_length=8

# list all paths that we want to compute:
srun generate_all_path_combinations_for_given_length.sh "$max_path_length"

# generate the one step anueploidy matrix that will be used to compute path likelihoods:
srun create_one_step_anueploidy_markov_matrix.sh "$max_path_length"

# generate the genome doubling matrix:
srun create_one_step_genome_doubling_markov_matrix.sh "$max_path_length"

# using srun here means that we will wait for everything to finish before moving on

max_CN=$((2 ** $max_path_length))
path_description="p"$max_path_length"_v4"

max_jobs=400

for rep in {0..0}
do
  for p_up in {0..100}
  do
    for p_down in {0..100}
    do
      if [ $(($p_up+$p_down)) -lt 101 ]
      then
        FILE="MATRICES/subbed_u"$p_up"_d"$p_down"_"$path_description".csv"
        if ! [ -f $FILE ]
        then
          echo $FILE
          echo "doesn't exist!"
          
          num_jobs=$(squeue -u "$USER" | awk 'NR > 1' | wc -l)
          while [ $num_jobs -ge $max_jobs ]
          do
            sleep 10
            num_jobs=$(squeue -u "$USER" | awk 'NR > 1' | wc -l)
          done
          
          sbatch ./generate_path_likelihoods.sh $p_up $p_down $max_CN $max_path_length $path_description
        fi
      fi
      echo $p_up
      echo $p_down
    done
  done
done

python check_all_files_exist.py


