#!/bin/bash

#SBATCH --job-name=meta_analysis
#SBATCH --output=meta_analysis_%j.out
#SBATCH --error=meta_analysis_%j.err
#SBATCH --time=24:00:00
#SBATCH --partition=regular
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=4G

#module load python/3.9

python do_meta_analysis.py standard_simulation

