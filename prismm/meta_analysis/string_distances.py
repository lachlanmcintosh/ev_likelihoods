from typing import List, Optional, Dict, Union

def levenshtein_distance(first_str: str, second_str: str) -> int:
    """
    Calculates the Levenshtein distance between the first string and the second string.

    Args:
        first_str: The first string.
        second_str: The second string.

    Returns:
        The Levenshtein distance.
    """
    # If the strings are identical, return 0
    if first_str == second_str:
        return 0

    # Swap if the second string is shorter than the first
    if len(first_str) < len(second_str):
        first_str, second_str = second_str, first_str

    # If the first string is empty, return the length of the second string
    if not first_str:
        return len(second_str)

    previous_row = range(len(second_str) + 1)
    for i, column1 in enumerate(first_str):
        current_row = [i + 1]
        for j, column2 in enumerate(second_str):
            insertions = previous_row[j + 1] + 1
            deletions = current_row[j] + 1
            substitutions = previous_row[j] + (column1 != column2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row

    return previous_row[-1]


def calculate_edit_distances(results: List[Optional[Dict[str, Union[str, int]]]]) -> List[int]:
    """
    Calculates the edit distance for the given results.

    Args:
        results: The result set for which the edit distance is calculated.

    Returns:
        A list of calculated edit distances.
    """
    edit_distances = []
    for result in results:
        if result is None:
            continue
        genome_doublings = result['genome_doublings']
        est_genome_doublings = result['est_genome_doublings']
        edit_distance = levenshtein_distance(str(genome_doublings), str(est_genome_doublings))
        edit_distances.append(edit_distance)

    return edit_distances


def print_different_ev_strings(results: List[Optional[Dict[str, str]]]) -> None:
    """
    Prints the different Evolutionary Value strings from the results.

    Args:
        results: The results set to check for different Evolutionary Value strings.

    Returns:
        None
    """
    for result in results:
        if result is None:
            continue
        ev_string = result['ev_str']
        est_ev_string = result['est_ev_str']

        if ev_string != est_ev_string:
            print(f"ev_string: {ev_string}, est_ev_string: {est_ev_string}")

