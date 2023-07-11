import copy
from typing import Dict

def create_new_chromosome(old_chromosome: Dict, chrom_count: int, epoch: int) -> Dict:
    """
    Create a new chromosome by deep copying an old chromosome and updating
    its unique identifier, epoch created, and parent details.

    Args:
        old_chromosome (Dict): The chromosome to be copied.
        chrom_count (int): The current chromosome count.
        epoch (int): The current epoch.

    Returns:
        new_chromosome (Dict): The newly created chromosome.
    """
    assert isinstance(old_chromosome, dict), f"Expected dict, got {type(old_chromosome)}"
    assert isinstance(chrom_count, int), f"Expected int, got {type(chrom_count)}"
    assert isinstance(epoch, int), f"Expected int, got {type(epoch)}"

    new_chromosome = copy.deepcopy(old_chromosome)
    chrom_count += 1
    new_chromosome["unique_identifier"] = chrom_count
    new_chromosome["epoch_created"] = epoch
    new_chromosome["parent"] = old_chromosome["unique_identifier"]

    return new_chromosome
