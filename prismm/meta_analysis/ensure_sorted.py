import random
from typing import List, Dict, Tuple, Any

def is_sorted_by_worst_aic(solutions: List[Dict[str, Any]]) -> bool:
    """
    Check if the given list of solutions is sorted in descending order by their 'AIC' value.
    
    Args:
        solutions: A list of solutions. Each solution is a dictionary that includes 'AIC' key.

    Returns:
        True if the solutions list is sorted by 'AIC' in descending order, False otherwise.
    """

    # If the list is empty, it's trivially sorted
    if not solutions:
        return True

    # AIC values should decrease or remain the same as we move through the list
    for i in range(len(solutions) - 1):
        if solutions[i]['AIC'] < solutions[i + 1]['AIC']:
            # We found a pair out of order
            return False

    # We didn't find any pairs out of order
    return True

# An example use-case of the function
if __name__ == "__main__":
    solutions = [
        {"AIC": 10, "solution": "sol1"},
        {"AIC": 8, "solution": "sol2"},
        {"AIC": 6, "solution": "sol3"}
    ]

    print(is_sorted_by_worst_aic(solutions))  # Expected Output: True

