#!/bin/bash

echo "Directory contents:"
ls -alhrt /vast/scratch/users/lmcintosh/CN_SV_SNV_evolutionary_likelihood/SIMULATIONS/ | awk '{print "/vast/scratch/users/lmcintosh/CN_SV_SNV_evolutionary_likelihood/SIMULATIONS/"$NF}'
ls -alhrt /vast/scratch/users/lmcintosh/CN_SV_SNV_evolutionary_likelihood/SIMULATIONS/

log_files="/vast/scratch/users/lmcintosh/CN_SV_SNV_evolutionary_likelihood/SIMULATIONS/simulation_slurm-*.log"
for file in $log_files
do
    echo "File: $file"
    echo "Tail:"
    tail -n 50 "$file"
    echo " "
done

