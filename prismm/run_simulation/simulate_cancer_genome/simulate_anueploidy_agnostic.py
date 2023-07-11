from typing import List, Dict, Tuple
import random

from clonal_trees.run_simulation.simulate_cancer_genome.create_new_chromosome import create_new_chromosome

def simulate_agnostic_chromosome_loss(chromosomes: List[Dict], chrom_count: int) -> int:
    """
    Simulate the loss of a random number of chromosomes.

    :param chromosomes: A list of chromosomes represented as dictionaries.
    :param chrom_count: The current chromosome count.
    :return: The updated chromosome count after the loss.
    """
    assert isinstance(chromosomes, list), f"Expected list, got {type(chromosomes)}"
    assert isinstance(chrom_count, int), f"Expected int, got {type(chrom_count)}"

    for which in random.sample(range(len(chromosomes)), random.randint(0, len(chromosomes))):
        chromosomes[which]["dead"] = True
        chrom_count -= 1

    return chrom_count

def simulate_agnostic_chromosome_gain(chromosomes: List[Dict], epoch: int, chrom_count: int) -> int:
    """
    Simulate the gain of a random number of chromosomes.

    :param chromosomes: A list of chromosomes represented as dictionaries.
    :param epoch: The current epoch.
    :param chrom_count: The current chromosome count.
    :return: The updated chromosome count after the gain.
    """
    assert isinstance(chromosomes, list), f"Expected list, got {type(chromosomes)}"
    assert isinstance(epoch, int), f"Expected int, got {type(epoch)}"
    assert isinstance(chrom_count, int), f"Expected int, got {type(chrom_count)}"

    for which in random.choices(range(len(chromosomes)), k=random.randint(0, len(chromosomes))):
        if not chromosomes[which]["dead"]:
            new_chrom = create_new_chromosome(chromosomes[which], chrom_count, epoch)
            chromosomes.append(new_chrom)
            chrom_count += 1

    return chrom_count

def simulate_anueploidy_agnostic(
    simulated_chromosomes: Dict[str, List[Dict]],
    epoch: int,
    chrom_count: int
) -> Tuple[int, Dict[str, List[Dict]]]:
    """
    Simulate Anueploidy on a given set of chromosomes, agnostic of chromosome details.

    Randomly select chromosomes to lose and gain, updating their status accordingly.

    :param simulated_chromosomes: A dictionary with chromosome types as keys and a list of chromosomes as values.
                                  Each chromosome is represented as a dictionary.
    :param epoch: The current epoch.
    :param chrom_count: The current chromosome count.
    :return: The updated chromosome count and the updated chromosomes after simulation.
    """
    assert isinstance(simulated_chromosomes, dict), f"Expected dict, got {type(simulated_chromosomes)}"
    assert isinstance(epoch, int), f"Expected int, got {type(epoch)}"
    assert isinstance(chrom_count, int), f"Expected int, got {type(chrom_count)}"

    for chrom_type, chromosomes in simulated_chromosomes.items():
        chrom_count = simulate_agnostic_chromosome_loss(chromosomes, chrom_count)
        chrom_count = simulate_agnostic_chromosome_gain(chromosomes, epoch, chrom_count)

    return chrom_count, simulated_chromosomes
