"""
This module contains functions for checking various aspects of the simulated chromosomes.

The checks include counting the number of chromosomes of a particular chromosome type 
with a specific paternal type that are not marked as dead, checking if all chromosomes 
in the simulated_chromosomes dictionary have unique identifiers, checking if all keys 
that are expected to be present in simulated chromosomes are actually present, and 
conducting several checks on the simulated chromosomes.
"""

from typing import List, Dict, Any, Union

def count_paternity(chromosomes: List[Dict[str, Any]], paternal: Union[str, bool]) -> int:
    """
    Count the number of chromosomes of a particular chromosome type with a specific paternal type that are not marked as dead.
    """

    for chrom in chromosomes:
        if "paternal" not in chrom:
            raise ValueError("Chromosome is missing 'paternal' key")
        if "dead" not in chrom:
            raise ValueError("Chromosome is missing 'dead' key")

    return sum(1 for chrom in chromosomes if chrom["paternal"] == paternal and not chrom["dead"])


def check_all_chrs_are_unique(simulated_chromosomes: Dict[str, List[Dict[str, str]]]) -> bool:
    """
    Check if all chromosomes in the simulated_chromosomes dictionary have unique identifiers.
    """

    ids = []
    for chrom_type in simulated_chromosomes.values():
        for chrom in chrom_type:
            if "unique_identifier" not in chrom:
                raise KeyError(f"The key 'unique_identifier' does not exist in {chrom}")
            ids.append(chrom["unique_identifier"])

    unique_ids = set(ids)
    if len(ids) != len(unique_ids):
        failed_ids = [id_ for id_ in ids if ids.count(id_) > 1]
        raise ValueError(f"Non-unique identifiers found: {failed_ids}")

    return True


def check_expected_keys_in_simulated_chromosomes_present(simulated_chromosomes: Dict[str, List[Dict[str, str]]]) -> bool:
    """
    Check if all keys that are expected to be present in simulated chromosomes are actually present.
    """

    expected_keys = {"SNVs", "paternal", "epoch_created", "parent", "unique_identifier", "dead"}

    for chrom_type in simulated_chromosomes.values():
        for chrom in chrom_type:
            if not expected_keys.issubset(chrom.keys()):
                return False

    return True


def check_simulated_chromosomes(simulated_chromosomes, pre, mid, post, ev_sequence):
    assert(pre + mid + post + 2 == len(ev_sequence))
    # some quick final sanity checks:
    if post == 0 or (mid == 0 and post == -1):
        assert(ev_sequence[-1] == "G")
        # if a genome doubling round just occurred then the copy number of every chromosome will be even:
        for chrom_type in simulated_chromosomes:
            paternity_count = count_paternity(simulated_chromosomes[chrom_type], paternal=True)
            assert paternity_count % 2 == 0, f"Chromosomes: {simulated_chromosomes[chrom_type]}, Paternity count: {paternity_count}"

            paternity_count = count_paternity(simulated_chromosomes[chrom_type], paternal=False)
            assert paternity_count % 2 == 0, f"Chromosomes: {simulated_chromosomes[chrom_type]}, Paternity count: {paternity_count}"

    for chrom_type in simulated_chromosomes:
        assert(len(simulated_chromosomes[chrom_type]) != 0)

    # some basic checks about the simulated chromosomes:
    assert(check_all_chrs_are_unique(simulated_chromosomes))
    assert(check_expected_keys_in_simulated_chromosomes_present(simulated_chromosomes))
