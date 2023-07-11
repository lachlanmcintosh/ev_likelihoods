from typing import Dict, List
from collections import Counter

def count_copy_numbers(simulated_chromosomes: Dict[str, List[Dict]]) -> Dict[str, List[int]]:
    """
    Count the number of each parental specific copy number found in the genome.

    Args:
        simulated_chromosomes: The input simulated chromosomes.

    Returns:
        A dictionary mapping chromosome type to a list of counts of living chromosomes for each parent.
    """
    observed_copy_numbers = {}
    for chrom_type in simulated_chromosomes:
        observed_copy_numbers[chrom_type] = [
            len([x for x in simulated_chromosomes[chrom_type] if paternal == x["paternal"] and not x["dead"]])
            for paternal in [True, False]
        ]
    return observed_copy_numbers


def count_CN_multiplicities(observed_CNs: Dict[str, List[int]]) -> Dict[int, int]:
    """
    Count the multiplicities of each observed copy number in the genome.

    Args:
        observed_CNs: The input observed copy numbers.

    Returns:
        A dictionary mapping each observed copy number to its multiplicity.
    """
    multiplicities = Counter()

    for CN in observed_CNs.values():
        multiplicities.update(CN)

    return multiplicities


def count_copy_number_multiplicities(observed_copy_numbers: Dict[str, List[int]]) -> Dict[int, int]:
    """
    Count the multiplicities of each observed copy number in the genome.

    Args:
        observed_copy_numbers: The input observed copy numbers.

    Returns:
        A dictionary mapping each observed copy number to its multiplicity.
    """
    multiplicities = count_CN_multiplicities(observed_copy_numbers)

    return multiplicities
