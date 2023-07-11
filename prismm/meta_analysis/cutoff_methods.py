import logging
import numpy as np
from sklearn.metrics import accuracy_score
from statistics_operations import print_summary_statistics
from typing import Dict, List, Tuple, Optional
from collections import Counter

# Setting up logging
logging.basicConfig(level=logging.INFO)

def calculate_average_tcn(best_estimate: Dict) -> float:
    """
    Calculates the average total copy number (TCN) for a given best estimate.
    
    Args:
        best_estimate (dict): The best estimate dictionary.

    Returns:
        float: The average TCN.
    """
    assert isinstance(best_estimate, dict), 'best_estimate should be a dictionary.'
    
    tcn = [best_estimate[chrom]["tree"][1][0] + best_estimate[chrom]["tree"][2][0] for chrom in range(23)]
    return np.mean(tcn)

def calculate_tcn_per_estimate(best_estimates: List[Dict]) -> List[float]:
    """
    Calculates the average TCN for each best estimate in the given list.
    
    Args:
        best_estimates (list): List of best estimate dictionaries.

    Returns:
        list: List of average TCNs for each best estimate.
    """
    assert isinstance(best_estimates, list), 'best_estimates should be a list.'
    
    return [calculate_average_tcn(best_estimate) for best_estimate in best_estimates]

def print_statistics_for_genome_doubling(best_estimates: List[Dict], TCNs: List[float], genome_doubling: int) -> None:
    """
    Prints the summary statistics for the given number of genome doublings.
    
    Args:
        best_estimates (list): List of best estimate dictionaries.
        TCNs (list): List of average TCNs for each best estimate.
        genome_doubling (int): The number of genome doublings.

    Returns:
        None
    """
    assert isinstance(best_estimates, list), 'best_estimates should be a list.'
    assert isinstance(TCNs, list), 'TCNs should be a list.'
    assert isinstance(genome_doubling, int), 'genome_doubling should be an integer.'
    
    filtered_tcn = [TCN for best_estimate, TCN in zip(best_estimates, TCNs) if best_estimate['genome_doublings'] == genome_doubling]
    data_pairs = list(zip(filtered_tcn, [genome_doubling]*len(filtered_tcn)))
    logging.info(f"Average TCN for estimates with {genome_doubling} genome doublings:")
    summary = print_summary_statistics(data_pairs,'TCN')
    logging.info(summary)
    return summary



def find_best_cutoff(best_estimates: List[Dict], num_cutoffs: int = 100) -> Tuple[np.ndarray, np.ndarray, np.ndarray, Tuple[float, float]]:
    """
    Finds the best cutoff values for the given list of best estimates.
    
    Args:
        best_estimates (list): List of best estimate dictionaries.
        num_cutoffs (int, optional): Number of cutoffs to consider. Defaults to 100.

    Returns:
        tuple: Tuple containing two 1D numpy arrays representing the cutoff values, 
        a 2D numpy array representing the accuracy values, and a tuple of the best cutoff values.
    """
    assert isinstance(best_estimates, list), 'best_estimates should be a list.'
    assert isinstance(num_cutoffs, int), 'num_cutoffs should be an integer.'
    
    TCNs = calculate_tcn_per_estimate(best_estimates)

    tcn_summaries = {}
    for gd in range(3):
        summary = print_statistics_for_genome_doubling(best_estimates, TCNs, gd)
        tcn_summaries[f"gd_{gd}_avg"] = summary

    min_TCN, max_TCN = min(TCNs), max(TCNs)
    cutoff_values = np.linspace(min_TCN, max_TCN, num_cutoffs)

    best_accuracy, best_cutoff_1, best_cutoff_2 = 0, None, None
    accuracy_values = np.full((num_cutoffs, num_cutoffs), np.nan)

    actual_gd = [best_estimate['genome_doublings'] for best_estimate in best_estimates]

    for i, cutoff_1 in enumerate(cutoff_values):
        for j in range(i+1, len(cutoff_values)):
            cutoff_2 = cutoff_values[j]

            predicted_gd = [
                2 if TCN > cutoff_2 else
                1 if TCN > cutoff_1 else
                0
                for TCN in TCNs
            ]

            accuracy = accuracy_score(actual_gd, predicted_gd)

            if accuracy > best_accuracy:
                best_accuracy, best_cutoff_1, best_cutoff_2 = accuracy, cutoff_1, cutoff_2

            accuracy_values[i, j] = accuracy

    logging.info(f"Best cutoffs: {best_cutoff_1}, {best_cutoff_2}, Accuracy: {best_accuracy}")

    return cutoff_values, cutoff_values, accuracy_values, (best_cutoff_1, best_cutoff_2), tcn_summaries

def find_best_accuracy_box(cutoff1_values: np.ndarray, cutoff2_values: np.ndarray, accuracy_values: np.ndarray) -> Tuple[float, float, float, float, float]:
    """
    Finds the accuracy box with the best accuracy.
    
    Args:
        cutoff1_values (np.ndarray): 1D numpy array representing the cutoff1 values.
        cutoff2_values (np.ndarray): 1D numpy array representing the cutoff2 values.
        accuracy_values (np.ndarray): 2D numpy array representing the accuracy values.

    Returns:
        tuple: Tuple of the minimum and maximum cutoff1 values, minimum and maximum cutoff2 values, and the best accuracy.
    """
    assert isinstance(cutoff1_values, np.ndarray), 'cutoff1_values should be a numpy array.'
    assert isinstance(cutoff2_values, np.ndarray), 'cutoff2_values should be a numpy array.'
    assert isinstance(accuracy_values, np.ndarray), 'accuracy_values should be a numpy array.'
    
    best_accuracy = np.nanmax(accuracy_values)
    cutoff1_list, cutoff2_list = [], []

    for i in range(len(cutoff1_values)):
        for j in range(len(cutoff2_values)):
            if accuracy_values[i][j] == best_accuracy:
                cutoff1_list.append(cutoff1_values[i])
                cutoff2_list.append(cutoff2_values[j])

    return min(cutoff1_list), max(cutoff1_list), min(cutoff2_list), max(cutoff2_list), best_accuracy


def calculate_accuracy_for_given_cutoffs(best_estimates: List[Dict], cutoff1: float, cutoff2: float) -> float:
    """
    Calculates the accuracy of the given cutoffs for the given best estimates.

    Args:
        best_estimates (list): List of best estimate dictionaries.
        cutoff1 (float): The first cutoff value.
        cutoff2 (float): The second cutoff value.

    Returns:
        float: The accuracy of the given cutoffs.
    """
    assert isinstance(best_estimates, list), 'best_estimates should be a list.'
    assert isinstance(cutoff1, (int, float)), 'cutoff1 should be a number.'
    assert isinstance(cutoff2, (int, float)), 'cutoff2 should be a number.'

    print(f"cutoff1: {cutoff1}, cutoff2: {cutoff2}")

    TCNs = calculate_tcn_per_estimate(best_estimates)
    actual_gd = [best_estimate['genome_doublings'] for best_estimate in best_estimates]

    predicted_gd = [
        2 if TCN > cutoff2 else
        1 if TCN > cutoff1 else
        0
        for TCN in TCNs
    ]

    print(f"Predicted gd counts: {predicted_gd}")
    print(f"Actual gd counts: {actual_gd}")

    accuracy = accuracy_score(actual_gd, predicted_gd)

    return accuracy

