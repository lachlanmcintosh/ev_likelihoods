import argparse
import pickle
import numpy as np
import glob
import os
import shelve
import sys
from general_functions import *
import re
import copy
import pandas as pd
from scipy import stats
from sklearn.metrics import accuracy_score
from typing import Any, Dict, List, Optional, Tuple, Union


def get_filenames_sorted_by_time(simulation_filename: str) -> List[str]:
    """
    Fetch all filenames associated with a given simulation and return them sorted by their creation times.

    Args:
        simulation_filename: The name of the simulation.

    Returns:
        A list of filenames sorted by their creation times.
    """
    file_pattern = f"SIMULATIONS/{simulation_filename}_*smaller.pickle"
    files = glob.glob(file_pattern)

    # Fetch creation times for all files
    creation_times = [os.path.getctime(file) for file in files]

    # Pair filenames with their creation times
    name_time_pairs = list(zip(files, creation_times))

    # Sort pairs by creation time
    name_time_pairs.sort(key=lambda x: x[1])

    # Unpack sorted filenames
    sorted_filenames = [pair[0] for pair in name_time_pairs]

    print(sorted_filenames)

    return sorted_filenames

def create_3x3_grid(dicts_list):
    grid = np.zeros((3, 3))
    for d in dicts_list:
        if d == None:
            continue
        gd = d['genome_doublings']
        gd_est = d['est_genome_doublings']
        if gd < 3 and gd_est < 3:
            grid[gd, gd_est] += 1
    return grid

def get_correct_proportion(grid):
    diagonal_sum = sum(grid[i][i] for i in range(len(grid)))
    total_sum = sum(sum(row) for row in grid)
    if total_sum == 0:
        return 0
    else:
        return diagonal_sum / total_sum

def get_all_best_estimates(list_of_lists):
    return [x[-1] for x in list_of_lists]
    
def load_file(filename: str, can_expect_EOF_error: bool = False):
    try:
        with open(filename, 'rb') as f:
            SS = pickle.load(f)
    except EOFError:
        if can_expect_EOF_error:
            print(f"EOFError: {filename} is incomplete. The file might still be generating.")
            return None
        else:
            raise
    except Exception as e:
        print(f"An error occurred while loading {filename}: {e}")
        return None

    return SS

def process_file(filename, can_expect_EOF_error):
    SS = load_file(filename, can_expect_EOF_error)

    if SS is None or "solutions" not in SS:
        return []

    print(filename)
    print("length of results:" + str(len(SS["solutions"])))
    for result in SS["solutions"]:
        test_case = int(re.search(r'(\d+)', filename).group())
        result['test_case'] = test_case

    for result in SS["solutions"]:
        print_result_info(result)

    return copy.deepcopy(SS["solutions"])


def process_files(filenames, can_expect_EOF_error):
    list_of_lists = []
    for filename in filenames:
        result = process_file(filename, can_expect_EOF_error)
        if result:
            ensure_sorted_by_worst_aic(result)
            list_of_lists.append(result)
    return list_of_lists


def print_result_info(result):
    fields = ["pre", "mid", "post", "p_up", "p_down", "plambda"]
    est_fields = ["est_"+ f for f in fields]

    print("(truth,estimated)")
    for field, est_field in zip(fields, est_fields):
        print(f"{field}: {result[field], result[est_field]}")

    print(f"AIC: {result['AIC']:.2f}")
    print("ev_string:" + str(result["ev_str"]))
    print("est_ev_string:" + str(result["est_ev_str"]))

import random
from typing import List, Dict, Tuple

def ensure_sorted_by_worst_aic(SS):
    if not SS["solutions"]:
        return True  # If the list is empty, it's trivially sorted
    for i in range(len(SS["solutions"]) - 1):
        # AIC values should decrease or remain the same as we move through the list
        if SS["solutions"][i]['AIC'] < SS["solutions"][i + 1]['AIC']:
            return False  # We found a pair out of order
    return True  # We didn't find any pairs out of order

def levenshtein_distance(a, b):
    """
    Calculates the Levenshtein distance between a and b.
    """
    size_x = len(a) + 1
    size_y = len(b) + 1
    matrix = [[0 for _ in range(size_y)] for _ in range(size_x)]
    for x in range(size_x):
        matrix [x][0] = x
    for y in range(size_y):
        matrix [0][y] = y

    for x in range(1, size_x):
        for y in range(1, size_y):
            if a[x-1] == b[y-1]:
                matrix [x][y] = min(
                    matrix[x-1][y] + 1,
                    matrix[x-1][y-1],
                    matrix[x][y-1] + 1
                )
            else:
                matrix [x][y] = min(
                    matrix[x-1][y] + 1,
                    matrix[x-1][y-1] + 1,
                    matrix[x][y-1] + 1
                )
    return matrix[size_x - 1][size_y - 1]


import collections
from statistics import mean, median, stdev

def calculate_edit_distances(results):
    distances = []
    for result in results:
        if result:
            ev_string = result['ev_str']
            est_ev_string = result['est_ev_str']
            distance = levenshtein_distance(ev_string, est_ev_string)
            distances.append(distance)
    return distances

def print_different_ev_strings(results):
    for result in results:
        if result:
            ev_string = result['ev_str']
            est_ev_string = result['est_ev_str']
            if ev_string != est_ev_string:
                print(f"ev_string: {ev_string}, est_ev_string: {est_ev_string}")


import numpy as np
import collections


def empirical_distribution(distances):
    rounded_distances = [round(d, 2) for d in distances]
    counter = collections.Counter(rounded_distances)
    distribution = {key: (val, round(val / len(distances), 3)) for key, val in counter.items()}
    distribution = {key: distribution[key] for key in sorted(distribution)}
    return distribution


import scipy.stats as stats

def print_summary_statistics(data_pairs, title):
    diffs = [abs(actual - estimated) for actual, estimated in data_pairs]
    diffs = np.array(diffs)
    
    mean = np.mean(diffs)
    std_err = stats.sem(diffs)
    confidence_interval = stats.t.interval(0.95, len(diffs)-1, loc=mean, scale=std_err)

    print(f"=== {title} Summary Statistics ===")
    print("\tMin:", np.min(diffs))
    print("\tMax:", np.max(diffs))
    print("\tMean:", mean)
    print("\tMedian:", np.median(diffs))
    print("\tMode:", stats.mode(diffs, keepdims=False)[0])
    print("\tStd Dev:", np.std(diffs))
    print("\t1st Quartile (25%):", np.percentile(diffs, 25))
    print("\t3rd Quartile (75%):", np.percentile(diffs, 75))
    print("\tInterquartile Range (IQR):", stats.iqr(diffs))
    print("\t95% Confidence Interval for Mean:", confidence_interval)
    print("\tProportion with zero difference:", np.mean(diffs == 0))

    print(f"\n=== {title} Empirical Distributions ===")
    dist = empirical_distribution(diffs)
    for key in sorted(dist):
        count, proportion = dist[key]
        print(f"\tDistance {key}: Proportion = {proportion}, Count = {count}")

    print("\n")



def summarize_p_similarity(all_results):
    p_up_pairs = [(result["p_up"], result["est_p_up"]/100.0) for result in all_results]
    p_down_pairs = [(result["p_down"], result["est_p_down"]/100.0) for result in all_results]

    print_summary_statistics(p_up_pairs, "p_up diffs")
    print_summary_statistics(p_down_pairs, "p_down diffs")

def summarize_p_similarity_zero_edit_distance(all_results):
    p_up_pairs = []
    p_down_pairs = []

    for result in all_results:
        edit_distance = levenshtein_distance(result['ev_str'], result['est_ev_str'])
        if edit_distance == 0:
            p_up_pairs.append((result["p_up"], result["est_p_up"]/100.0))
            p_down_pairs.append((result["p_down"], result["est_p_down"]/100.0))

    print_summary_statistics(p_up_pairs, "p_up diffs for zero edit distance")
    print_summary_statistics(p_down_pairs, "p_down diffs for zero edit distance")

def summarize_plambda_similarity_zero_edit_distance(all_results):
    plambda_pairs = []

    for result in all_results:
        if 'A' not in result['ev_str']:
            plambda_pairs.append((1, 1))  # zero difference
        else:
            edit_distance = levenshtein_distance(result['ev_str'], result['est_ev_str'])
            if edit_distance == 0:
                plambda_pairs.append((1, result["est_plambda"]/result["plambda"]))

    print_summary_statistics(plambda_pairs, "plambda diffs for zero edit distance")

def get_plambda_data_pairs(results):
    plambda_data_pairs = [(1, (result['est_plambda']/result["plambda"]) if 'A' in result['ev_str'] else 1) for result in results]
    return plambda_data_pairs

pd.set_option('display.width', 120)

def aggregate_similarity_scores(list_of_dfs):
    similarity_scores = [[] for _ in range(len(list_of_dfs[0]))]

    for df in list_of_dfs:
        for i, score in enumerate(df['Similarity Score']):
            # Only append the score if it is a number
            if pd.notnull(score):
                similarity_scores[i].append(score)

    min_values = [min(scores) if scores else np.nan for scores in similarity_scores]
    max_values = [max(scores) if scores else np.nan for scores in similarity_scores]
    mean_values = [np.mean(scores) if scores else np.nan for scores in similarity_scores]
    median_values = [np.median(scores) if scores else np.nan for scores in similarity_scores]
    std_dev_values = [np.std(scores, ddof=0) if scores else np.nan for scores in similarity_scores]
    
    mode_values = []
    for scores in similarity_scores:
        if scores and len(set(scores)) > 1:
            mode_score = stats.mode(scores, nan_policy='omit')[0]
            if np.isscalar(mode_score):
                mode_values.append(mode_score)
            else:
                mode_values.append(mode_score[0])
        else:
            mode_values.append(scores[0] if scores else np.nan)

    iqr_values = [stats.iqr(scores, nan_policy='omit') if scores else np.nan for scores in similarity_scores]

    agg_df = list_of_dfs[0].copy()
    agg_df = agg_df.drop('Similarity Score', axis=1)
    agg_df.rename(columns={
        'Copy Number Distance Function': 'CN dist',
        'Epoch Created Distance Function': 'EC dist',
        'Normalising Constant': 'Norm Const'}, inplace=True)
    agg_df['Min'] = min_values
    agg_df['Max'] = max_values
    agg_df['Mean'] = mean_values
    agg_df['Median'] = median_values
    agg_df['StdDev'] = std_dev_values
    agg_df['Mode'] = mode_values
    agg_df['IQR'] = iqr_values

    # Round values to 3 decimal places
    agg_df = agg_df.round(3)

    return agg_df


def calculate_TCNs(best_estimates):
    TCNs = []
    for best_estimate in best_estimates:
        tcn = []
        for chrom in range(23):
            tree = best_estimate[chrom]["tree"]
            tcn.append(tree[1][0] + tree[2][0])
        average_tcn = np.mean(tcn)
        TCNs.append(average_tcn)
    return TCNs



def find_best_cutoff(best_estimates, num_cutoffs=100):
    TCNs = calculate_TCNs(best_estimates)
   
    # Split TCNs based on 'genome_doublings'
    TCNs_0 = [TCN for best_estimate, TCN in zip(best_estimates, TCNs) if best_estimate['genome_doublings'] == 0]
    TCNs_1 = [TCN for best_estimate, TCN in zip(best_estimates, TCNs) if best_estimate['genome_doublings'] == 1]
    TCNs_2 = [TCN for best_estimate, TCN in zip(best_estimates, TCNs) if best_estimate['genome_doublings'] == 2]

    # Print summary statistics
    zero_vector = [0]*len(TCNs_0)
    data_pairs_0 = list(zip(TCNs_0, zero_vector))
    print_summary_statistics(data_pairs_0, "Average TCN for estimates with 0 genome doublings")

    one_vector = [1]*len(TCNs_1)
    data_pairs_1 = list(zip(TCNs_1, one_vector))
    print_summary_statistics(data_pairs_1, "Average TCN for estimates with 1 genome doubling")

    two_vector = [2]*len(TCNs_2)
    data_pairs_2 = list(zip(TCNs_2, two_vector))
    print_summary_statistics(data_pairs_2, "Average TCN for estimates with 2 genome doublings")


    min_TCN = min(TCNs)
    max_TCN = max(TCNs)

    cutoff_values = np.linspace(min_TCN, max_TCN, num_cutoffs)

    best_cutoff_1, best_cutoff_2 = None, None
    best_accuracy = 0

    # Initialize a 2D array for storing accuracies
    accuracy_values = np.empty((num_cutoffs, num_cutoffs))
    accuracy_values[:] = np.NaN

    # For each pair of cutoff values
    for i, cutoff_1 in enumerate(cutoff_values):
        for j, cutoff_2 in enumerate(cutoff_values):
            if cutoff_1 >= cutoff_2:
                continue  # we ensure cutoff_1 < cutoff_2

            # Predict GD category based on cutoffs
            predicted_gd = [
                2 if TCN > cutoff_2 else
                1 if TCN > cutoff_1 else
                0
                for TCN in TCNs
            ]

            # Actual GD category
            actual_gd = [
                best_estimate['genome_doublings']
                for best_estimate in best_estimates
            ]

            # Calculate accuracy
            accuracy = accuracy_score(actual_gd, predicted_gd)

            # Update best cutoffs if this pair gives a better accuracy
            if accuracy > best_accuracy:
                best_accuracy = accuracy
                best_cutoff_1 = cutoff_1
                best_cutoff_2 = cutoff_2

            # Store accuracy in matrix
            accuracy_values[i, j] = accuracy

    print(f"Best cutoffs: {best_cutoff_1}, {best_cutoff_2}, Accuracy: {best_accuracy}")

    return cutoff_values, cutoff_values, accuracy_values, (best_cutoff_1, best_cutoff_2)


def find_best_accuracy_box(cutoff1_values, cutoff2_values, accuracy_values):
    # Find the best accuracy
    best_accuracy = np.nanmax(accuracy_values)

    # Initialize lists to store the cutoff values that give the best accuracy
    cutoff1_list = []
    cutoff2_list = []

    # Loop through each element in the accuracy matrix
    for i in range(len(cutoff1_values)):
        for j in range(len(cutoff2_values)):
            # If the accuracy is the best accuracy, record the corresponding cutoff values
            if accuracy_values[i][j] == best_accuracy:
                cutoff1_list.append(cutoff1_values[i])
                cutoff2_list.append(cutoff2_values[j])

    # Calculate the min and max of the cutoff values
    cutoff1_min = min(cutoff1_list)
    cutoff1_max = max(cutoff1_list)
    cutoff2_min = min(cutoff2_list)
    cutoff2_max = max(cutoff2_list)

    # Return the box boundaries and the best accuracy
    return cutoff1_min, cutoff1_max, cutoff2_min, cutoff2_max, best_accuracy

def write_results_to_file(filename, cutoff1_min, cutoff1_max, cutoff2_min, cutoff2_max, best_accuracy):
    result = {
        "Best accuracy": best_accuracy,
        "Best accuracy box boundaries": {
            "Cutoff1_min": cutoff1_min,
            "Cutoff1_max": cutoff1_max,
            "Cutoff2_min": cutoff2_min,
            "Cutoff2_max": cutoff2_max,
        },
    }

    with open(filename, 'w') as f:
        json.dump(result, f, indent=4)


def main(simulation_filename: str):
    """
    Main function of the script. Responsible for calling other helper functions and orchestrating
    the entire code flow.
    """
    sorted_filenames = get_filenames_sorted_by_time(simulation_filename)
    print(sorted_filenames)

    list_of_test_cases = process_files(sorted_filenames, can_expect_EOF_error = True)

    print("BEST ESTIMATED")
    best_estimates = get_all_best_estimates(list_of_test_cases)
    grid = create_3x3_grid(best_estimates)
    print(grid)
    correct_proportion = get_correct_proportion(grid)
    print(f"The proportion of times GD is correctly guessed is {correct_proportion}")

    print("Summarise levenshtein path distance:")
    lv_distances = calculate_edit_distances(best_estimates)
    distance_pairs = zip(lv_distances*100,[0]*len(lv_distances))
    print_summary_statistics(distance_pairs, "Levenshtein Path Distance")

    print("Print all the pairs of ev strings that are different:")
    print_different_ev_strings(best_estimates)

    print("\n\nSummarise the similarity in probability estimates")
    summarize_p_similarity(best_estimates)

    print("Summarise the similarity in probability estimates for when the path difference is zero")
    summarize_p_similarity_zero_edit_distance(best_estimates)

    print("Summarise the similarity in the SNV rates:")
    plambda_data_pairs = get_plambda_data_pairs(best_estimates)
    print_summary_statistics(plambda_data_pairs, "Normalised SNV rate differences")

    print("Summarise the similarity in the SNV rates for when the path difference is zero:")
    summarize_plambda_similarity_zero_edit_distance(best_estimates)

    # Filter best_estimates to only include those with Levenshtein distance 0
    best_estimates_zero_levenshtein = [best_estimate for best_estimate, distance in zip(best_estimates, lv_distances) if distance == 0]

    pd.set_option('display.max_columns', 200)

    # Pass the list of dataframes to the aggregate_similarity_scores function
    print("Summarise the tree similarity statistics")
    agg_df = aggregate_similarity_scores([result['distance_dataframe'] for result in best_estimates])
    print(agg_df)

    print("Summarise the tree similarity statistics when Levenshtein distance is zero")
    agg_df_zero_levenshtein = aggregate_similarity_scores([result['distance_dataframe'] for result in best_estimates_zero_levenshtein])
    print(agg_df_zero_levenshtein)

    print("\nFinding the best cutoff in average Total Copy Number (TCN) for predicting 'genome_doublings':")

    cutoff1_values, cutoff2_values, accuracy_values, best_cutoffs = find_best_cutoff(best_estimates, num_cutoffs=20)
    np.set_printoptions(precision=3, suppress=True, threshold=np.inf, linewidth=200)
    print("Cutoff 1 values:", cutoff1_values)
    print("Cutoff 2 values:", cutoff2_values)
    print("Accuracy values:", accuracy_values)
    print("Best cutoffs:", best_cutoffs)

    cutoff1_min, cutoff1_max, cutoff2_min, cutoff2_max, best_accuracy = find_best_accuracy_box(cutoff1_values, cutoff2_values, accuracy_values)
    print("Best accuracy:", best_accuracy)
    print("Best accuracy box boundaries: Cutoff1_min = {}, Cutoff1_max = {}, Cutoff2_min = {}, Cutoff2_max = {}".format(cutoff1_min, cutoff1_max, cutoff2_min, cutoff2_max))

    result_filename = f"{simulation_filename}_results.json"
    write_results_to_file(result_filename, cutoff1_min, cutoff1_max, cutoff2_min, cutoff2_max, best_accuracy)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--simulation_filename", type=str, required=True, help="Name of the simulation")
    args = parser.parse_args()
    
    main(args.simulation_filename)
