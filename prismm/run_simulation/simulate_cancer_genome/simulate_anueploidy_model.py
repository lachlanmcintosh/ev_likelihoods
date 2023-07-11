from typing import Dict, List, Tuple
import copy
import numpy as np

from clonal_trees.run_simulation.simulate_cancer_genome.create_new_chromosome import create_new_chromosome

def simulate_chromosome_changes(
    chromosomes: List[Dict],
    epoch: int,
    chrom_count: int,
    p_up: float,
    p_down: float
) -> Tuple[List[Dict], int]:
    """
    Simulate changes in a list of chromosomes based on the provided probabilities.

    Args:
        chromosomes (List[Dict]): A list of chromosomes represented as dictionaries.
        epoch (int): The current epoch.
        chrom_count (int): The current chromosome count.
        p_up (float): The probability for a chromosome to be duplicated.
        p_down (float): The probability for a chromosome to be lost.

    Returns:
        new_chromosomes (List[Dict]): The list of chromosomes after the simulation.
        chrom_count (int): The updated chromosome count after the simulation.
    """
    new_chromosomes = []

    for chrom in copy.deepcopy(chromosomes):
        # need to use deepcopy here so that we don't overwrite the "deadness" of current chromosomes
        change = np.random.choice([1, 0, -1], 1, p=[p_up, 1 - p_up - p_down, p_down])[0]

        if chrom["dead"] or change == 0:
            new_chromosomes += [chrom]
        elif change == 1:
            new_chromosome = create_new_chromosome(chrom, chrom_count, epoch)
            new_chromosomes += [chrom, new_chromosome]
            chrom_count += 1
        elif change == -1:
            chrom["dead"] = True
            new_chromosomes += [chrom]

    return new_chromosomes, chrom_count

def simulate_anueploidy(
    simulated_chromosomes: Dict[str, List[Dict]],
    epoch: int,
    chrom_count: int,
    p_up: float,
    p_down: float
) -> Tuple[int, Dict[str, List[Dict]]]:
    """
    Simulate Anueploidy on a given set of chromosomes, taking into account specific
    probabilities for a chromosome to be duplicated or lost.

    Args:
        simulated_chromosomes (Dict[str, List[Dict]]): A dictionary with chromosome types as keys
        and a list of chromosomes as values. Each chromosome is represented as a dictionary.
        epoch (int): The current epoch.
        chrom_count (int): The current chromosome count.
        p_up (float): The probability for a chromosome to be duplicated.
        p_down (float): The probability for a chromosome to be lost.

    Returns:
        chrom_count (int): The updated chromosome count.
        simulated_chromosomes (Dict[str, List[Dict]]): The updated chromosomes after simulation.
    """
    assert isinstance(simulated_chromosomes, dict), f"Expected dict, got {type(simulated_chromosomes)}"
    assert isinstance(epoch, int), f"Expected int, got {type(epoch)}"
    assert isinstance(chrom_count, int), f"Expected int, got {type(chrom_count)}"
    assert isinstance(p_up, float), f"Expected float, got {type(p_up)}"
    assert isinstance(p_down, float), f"Expected float, got {type(p_down)}"
    assert 0 <= p_up <= 1, "p_up should be between 0 and 1"
    assert 0 <= p_down <= 1, "p_down should be between 0 and 1"
    assert p_up + p_down <= 1, "The sum of p_up and p_down should be less than or equal to 1"

    for chrom_type, chromosomes in simulated_chromosomes.items():
        original_chrom_count = chrom_count  # store the original count before the loop
        new_chromosomes, temp_chrom_count = simulate_chromosome_changes(chromosomes, epoch, chrom_count, p_up, p_down)
        while not any(chrom["dead"] == False for chrom in new_chromosomes):
            new_chromosomes, temp_chrom_count = simulate_chromosome_changes(chromosomes, epoch, original_chrom_count, p_up, p_down)
        simulated_chromosomes[chrom_type] = new_chromosomes
        chrom_count = temp_chrom_count  # update the count only after a successful loop

    return chrom_count, simulated_chromosomes
