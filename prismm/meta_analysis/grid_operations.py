import numpy as np
from typing import List, Optional, Dict

def create_3x3_grid(dicts_list: List[Optional[Dict[str, int]]]) -> np.ndarray:
    """
    Create a 3x3 grid from a list of dictionaries.

    Parameters:
    dicts_list (List[Optional[Dict[str, int]]]): List of dictionaries. Each dictionary should contain 'genome_doublings'
                                                  and 'est_genome_doublings' keys. The dictionaries can also be None.

    Returns:
    np.ndarray: 3x3 numpy array with the count of dictionaries where 'genome_doublings' and 'est_genome_doublings'
                are less than 3.
    """
    # Initialize 3x3 grid
    grid = np.zeros((3, 3), dtype=int)

    # Iterate over dictionaries in the list
    for dict_item in dicts_list:
        # Skip iteration if dictionary is None
        if dict_item is None:
            continue

        # Check if necessary keys are present in the dictionary
        if 'genome_doublings' not in dict_item or 'est_genome_doublings' not in dict_item:
            raise KeyError("'genome_doublings' or 'est_genome_doublings' missing from dictionary.")

        genome_doublings = dict_item['genome_doublings']
        est_genome_doublings = dict_item['est_genome_doublings']

        # Update grid if conditions are met
        if genome_doublings < 3 and est_genome_doublings < 3:
            grid[genome_doublings, est_genome_doublings] += 1

    return grid


def get_correct_proportion(grid: np.ndarray) -> float:
    """
    Compute the correct proportion from a 2D grid.

    Parameters:
    grid (np.ndarray): 2D numpy array.

    Returns:
    float: The proportion of the sum of the diagonal elements to the sum of all elements in the grid.
           If total_sum is zero, return zero to avoid division by zero error.
    """
    # Ensure the input is a 2D numpy array
    if not isinstance(grid, np.ndarray) or len(grid.shape) != 2:
        raise ValueError("Input should be a 2D numpy array.")

    diagonal_sum = grid.trace()
    total_sum = grid.sum()

    return diagonal_sum / total_sum if total_sum != 0 else 0

