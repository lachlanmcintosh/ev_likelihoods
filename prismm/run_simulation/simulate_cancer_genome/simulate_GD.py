from typing import Dict, List, Tuple

from clonal_trees.run_simulation.simulate_cancer_genome.create_new_chromosome import create_new_chromosome

def simulate_gd(
    simulated_chromosomes: Dict[str, List[Dict]],
    epoch: int,
    chrom_count: int
) -> Tuple[int, Dict[str, List[Dict]]]:
    """
    Simulate Genome Doubling (GD) on a given set of chromosomes.

    For each chromosome, if it is not dead, create a new chromosome by copying
    it and updating the unique identifier, epoch created, and parent details.

    Args:
        simulated_chromosomes (Dict[str, List[Dict]]): A dictionary with chromosome types as keys
        and a list of chromosomes as values. Each chromosome is represented as a dictionary.

        epoch (int): The current epoch.

        chrom_count (int): The current chromosome count.

    Returns:
        chrom_count (int): The updated chromosome count.

        simulated_chromosomes (Dict[str, List[Dict]]): The updated chromosomes after simulation.
    """
    assert isinstance(simulated_chromosomes, dict), f"Expected dict, got {type(simulated_chromosomes)}"
    assert isinstance(epoch, int), f"Expected int, got {type(epoch)}"
    assert isinstance(chrom_count, int), f"Expected int, got {type(chrom_count)}"

    for chrom_type, chromosomes in simulated_chromosomes.items():
        new_chromosomes = []
        for chrom in chromosomes:
            if not chrom["dead"]:
                new_chromosomes.append(create_new_chromosome(chrom, chrom_count, epoch))
                chrom_count += 1
        simulated_chromosomes[chrom_type] += new_chromosomes

    return chrom_count, simulated_chromosomes
