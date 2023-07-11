#!/bin/bash

simulation_name="find_best_parameter_simulation"

# Clean log and text directories
function clean_directories() {
    local simulation_name=$1
    local log_directory="/vast/scratch/users/lmcintosh/CN_SV_SNV_evolutionary_likelihood/SIMULATIONS/${simulation_name}_*.log"
    local txt_directory="/vast/scratch/users/lmcintosh/CN_SV_SNV_evolutionary_likelihood/SIMULATIONS/${simulation_name}_*.txt*"

    rm ${log_directory}
    rm ${txt_directory}
}

# Generate random parameters for a single job
function generate_parameters() {
    local solutions=$((RANDOM % 100 + 1))
    local default_path_length=$((RANDOM % 11))
    local prob_dist_filter=${prob_dist_filters[$((RANDOM % ${#prob_dist_filters[@]}))]}
    local path_length_diff=$((RANDOM % 11))

    echo $solutions $default_path_length $prob_dist_filter $path_length_diff
}

# Submit a single job with the given parameters
function submit_single_job() {
    local i=$1
    local solutions=$2
    local default_path_length=$3
    local prob_dist_filter=$4
    local path_length_diff=$5
    local simulation_name=$6

    sbatch \
    --job-name=singlecore_job \
    --ntasks=1 \
    --mem=40G \
    --requeue \
    --qos=bonus \
    --output=/vast/scratch/users/lmcintosh/CN_SV_SNV_evolutionary_likelihood/SIMULATIONS/${simulation_name}_${i}_%j.log \
    --error=/vast/scratch/users/lmcintosh/CN_SV_SNV_evolutionary_likelihood/SIMULATIONS/${simulation_name}_${i}_%j.log \
    --wrap "python do_all.py $i -s $solutions -d $default_path_length -p $prob_dist_filter -l $path_length_diff -f $simulation_name --run_scripts 1 2"
}

# Submit all jobs with a delay between batches
function submit_jobs() {
    local atatime=$1
    local reps=$2

    for rep in $(seq 1 $reps)
    do
        for i in $(seq $((atatime * (rep - 1) + 1)) $((atatime * rep)))
        do
            local parameters=($(generate_parameters))
            submit_single_job $i ${parameters[@]} $simulation_name
        done

        # Sleep for 2 seconds between batches
        sleep 2
    done
}

# Run touch_all.sh
sbatch touch_all.sh

# Clean directories
clean_directories $simulation_name

# Set parameters for job submission
atatime=40
reps=12
prob_dist_filters=(1.00000000e-05 3.16227766e-05 1.00000000e-04 3.16227766e-04 1.00000000e-03 3.16227766e-03 1.00000000e-02 3.16227766e-02 1.00000000e-01 3.16227766e-01 1.00000000e+00)

# Submit the array job
submit_jobs $atatime $reps

