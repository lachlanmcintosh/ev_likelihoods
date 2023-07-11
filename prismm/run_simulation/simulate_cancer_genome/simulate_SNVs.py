from typing import Dict, List, Any, Tuple
import numpy as np

def simulate_snvs(simulated_chromosomes: Dict[int, List[Dict[str, Any]]],
                  lengths: Dict[int, float],
                  rate: float,
                  epoch: int,
                  snv_count: int) -> Tuple[int, Dict[int, List[Dict[str, Any]]]]:
    """
    Simulate SNVs (Single Nucleotide Variants) and update the simulated_chromosomes.

    Args:
        simulated_chromosomes: A dictionary containing lists of chromosomes.
        lengths: A list of chromosome lengths.
        rate: The rate parameter for the Poisson distribution.
        epoch: The current epoch.
        snv_count: The current count of SNVs.

    Returns:
        A tuple containing the updated SNV count and the updated simulated_chromosomes.
    """
    for chrom_type, chrom_list in simulated_chromosomes.items():
        for chrom in chrom_list:
            if chrom["dead"]:
                continue

            additional_snv_count = np.random.poisson(rate * lengths[chrom_type])

            for x in range(snv_count + 1, snv_count + additional_snv_count + 1):
                chrom["SNVs"].append({"unique_identifier": str(x), "epoch_created": epoch})

            snv_count += additional_snv_count

    return snv_count, simulated_chromosomes

    
