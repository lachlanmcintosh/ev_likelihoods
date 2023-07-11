from typing import Dict, List

def initialize_simulated_chromosomes() -> Dict[int, List[Dict[str, str]]]:
    """
    Initialize a dictionary of simulated chromosomes with default values.

    Returns:
        A dictionary containing lists of initialized chromosomes.
    """
    simulated_chromosomes = {}

    for chrom_type in range(23):
        simulated_chromosomes[chrom_type] = [
            {
                "unique_identifier": chrom_type + x,
                "parent": -1,  # We use -1 to denote the root of the tree.
                "epoch_created": 0,
                "paternal": bool(x == 0),  # Each chromosome is paternal or not.
                "SNVs": [],
                "dead": False,
            } for x in (0,23)
        ]

    return simulated_chromosomes
