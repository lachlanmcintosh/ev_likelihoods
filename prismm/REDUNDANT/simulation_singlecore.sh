#!/bin/bash
#SBATCH --job-name=singlecore_job
#SBATCH --ntasks=1
#SBATCH --mem=40G
#SBATCH --requeue
#SBATCH --qos=bonus
#SBATCH --array=1-<MAX_INDEX>
#SBATCH --output=/vast/scratch/users/lmcintosh/CN_SV_SNV_evolutionary_likelihood/SIMULATIONS/simulation_slurm-%j_%a.log
#SBATCH --error=/vast/scratch/users/lmcintosh/CN_SV_SNV_evolutionary_likelihood/SIMULATIONS/simulation_slurm-%j_%a.log


source ~/anaconda3/etc/profile.d/conda.sh
conda activate sage

#python run_simulation_and_analysis2.py $SLURM_ARRAY_TASK_ID
python do_all.py $SLURM_ARRAY_TASK_ID -f "GD_tree_simulation" --run_scripts 1 2 3 4 
