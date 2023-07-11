#!/bin/bash
#SBATCH --job-name=singlecore_job 
#SBATCH --ntasks=1
#SBATCH --mem=40G
#SBATCH --requeue
#SBATCH --qos=bonus

sage everything.sage ${1} ${2} ${3} ${4} ${5}

