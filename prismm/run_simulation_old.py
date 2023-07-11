##### WRITE A FUNCTION THAT CAN SIMULATE A GENOME
import math
from constants import *
import numpy as np
import random
import copy
import sys
import pickle
from general_functions import *
import argparse
import logging
import json
from typing import Dict, List, Any, Tuple, Optional
from collections import Counter


import subprocess
import os



def count_paternity(chromosomes: List[Dict[str, any]], paternal: str) -> int:
    """
    Count the number of chromosomes of a particular chromosome type with a specific paternal type that are not marked as dead.

    Args:
        chromosomes: A list of dictionaries representing chromosomes of a particular chromosomal type/number. Each chromosome dictionary 
            must have a 'paternal' key indicating the type of paternal chromosome and a 'dead' key 
            indicating whether the chromosome is marked as dead or not.
        paternal: The type of paternal chromosome to count.

    Raises:
        ValueError: If a chromosome is missing the 'paternal' key or the 'dead' key.

    Returns:
        The number of chromosomes of the specified type with the specified paternal type that are not marked as dead.
    """ 
    for chrom in chromosomes:
        if "paternal" not in chrom:
            raise ValueError("Chromosome is missing 'paternal' key")
        if "dead" not in chrom:
            raise ValueError("Chromosome is missing 'dead' key")
    
    return sum(1 for chrom in chromosomes if chrom["paternal"] == paternal and not chrom["dead"])


def test_count_paternity():
    Input1 = ([
        {
            'unique_identifier': 22,
            'parent': -1,
            'epoch_created': 0,
            'paternal': True,
            'SNVs': [],
            'dead': False
        },
        {
            'unique_identifier': 45,
            'parent': -1,
            'epoch_created': 0,
            'paternal': False,
            'SNVs': [],
            'dead': False
        },
        {
            'unique_identifier': 91,
            'parent': 22,
            'epoch_created': 0,
            'paternal': True,
            'SNVs': [],
            'dead': False
        },
        {
            'unique_identifier': 92,
            'parent': 45,
            'epoch_created': 0,
            'paternal': False,
            'SNVs': [],
            'dead': False
        }
    ], False)
    Output1 = 2

    Input2 = ([
        {
            'unique_identifier': 19,
            'parent': -1,
            'epoch_created': 0,
            'paternal': True,
            'SNVs': [{'unique_identifier': '31', 'epoch_created': 2}],
            'dead': False
        },
        {
            'unique_identifier': 42,
            'parent': -1,
            'epoch_created': 0,
            'paternal': False,
            'SNVs': [],
            'dead': True
        },
        {
            'unique_identifier': 85,
            'parent': 19,
            'epoch_created': 0,
            'paternal': True,
            'SNVs': [],
            'dead': False
        },
        {
            'unique_identifier': 106,
            'parent': 86,
            'epoch_created': 1,
            'paternal': False,
            'SNVs': [],
            'dead': False
        },
        {
            'unique_identifier': 86,
            'parent': 42,
            'epoch_created': 0,
            'paternal': False,
            'SNVs': [],
            'dead': False
        },
        {
            'unique_identifier': 189,
            'parent': 19,
            'epoch_created': 2,
            'paternal': True,
            'SNVs': [{'unique_identifier': '31', 'epoch_created': 2}],
            'dead': False
        },
        {
            'unique_identifier': 190,
            'parent': 85,
            'epoch_created': 2,
            'paternal': True,
            'SNVs': [],
            'dead': False
        },
        {
            'unique_identifier': 191,
            'parent': 106,
            'epoch_created': 2,
            'paternal': False,
            'SNVs': [],
            'dead': False
        },
        {
            'unique_identifier': 192,
            'parent': 86,
            'epoch_created': 2,
            'paternal': False,
            'SNVs': [],
            'dead': False
        }
    ], True)
    Output2 = 4

    assert count_paternity(Input1[0], Input1[1]) == Output1
    assert count_paternity(Input2[0], Input2[1]) == Output2


test_count_paternity()


def check_all_chrs_are_unique(simulated_chromosomes: Dict[str, List[Dict[str, str]]]) -> bool:
    """
    Check if all chromosomes in the simulated_chromosomes dictionary have unique identifiers.

    Args:
        simulated_chromosomes: A dictionary containing lists of chromosomes.

    Returns:
        True if all chromosomes have unique identifiers.

    Raises:
        KeyError: If a chromosome does not have the 'unique_identifier' key.
        ValueError: If non-unique identifiers are found in the chromosomes.
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

import pytest
def test_check_all_chrs_are_unique():
    """
    Test the check_all_chrs_are_unique function.
    """
    unique_chrs = {
        "type_A": [{"unique_identifier": 1}, {"unique_identifier": 2}],
        "type_B": [{"unique_identifier": 3}, {"unique_identifier": 4}],
    }
    assert check_all_chrs_are_unique(unique_chrs)

    non_unique_chrs = {
        "type_A": [{"unique_identifier": 1}, {"unique_identifier": 2}],
        "type_B": [{"unique_identifier": 2}, {"unique_identifier": 4}],
    }
    with pytest.raises(ValueError):
        check_all_chrs_are_unique(non_unique_chrs)

    # More complex test cases
    non_unique_chrs_complex = {
        "type_A": [{"unique_identifier": i} for i in range(10)],
        "type_B": [{"unique_identifier": i} for i in range(5, 15)],
    }
    with pytest.raises(ValueError):
        check_all_chrs_are_unique(non_unique_chrs_complex)

    non_unique_chrs_single_type = {
        "type_A": [{"unique_identifier": 1}, {"unique_identifier": 1}],
    }
    with pytest.raises(ValueError):
        check_all_chrs_are_unique(non_unique_chrs_single_type)


test_check_all_chrs_are_unique()


def check_expected_keys_in_simulated_chromosomes_present(simulated_chromosomes: Dict[str, List[Dict[str, str]]]) -> bool:
    """
    Check if all keys that are expected to be present in simulated chromosomes are actually present.

    Args:
        simulated_chromosomes: A dictionary containing lists of chromosomes.

    Returns:
        True if all expected keys are present in all chromosomes, False otherwise.
    """
    expected_keys = {"SNVs", "paternal", "epoch_created", "parent", "unique_identifier", "dead"}
    
    for chrom_type in simulated_chromosomes.values():
        for chrom in chrom_type:
            if not expected_keys.issubset(chrom.keys()):
                return False
    
    return True

def test_check_expected_keys_in_simulated_chromosomes_present():
    valid_chrs = {
        "type_A": [
            {
                "SNVs": [],
                "paternal": "A",
                "epoch_created": 0,
                "parent": None,
                "unique_identifier": 1,
                "dead": False,
            }
        ],
    }
    assert check_expected_keys_in_simulated_chromosomes_present(valid_chrs)

    invalid_chrs = {
        "type_A": [
            {
                "SNVs": [],
                "paternal": "A",
                "epoch_created": 0,
                "unique_identifier": 1,
            }
        ],
    }
    assert not check_expected_keys_in_simulated_chromosomes_present(invalid_chrs)

test_check_expected_keys_in_simulated_chromosomes_present()


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

def test_initialize_simulated_chromosomes():
    chromosomes = initialize_simulated_chromosomes()

    assert len(chromosomes) == 23, "Number of chromosomes should be 23"

    for chrom_type, chrom_data in chromosomes.items():
        assert len(chrom_data) == 2, "Each chromosome type should have 2 entries"

        for entry in chrom_data:
            assert entry["unique_identifier"] == chrom_type + 23*(1-int(entry["paternal"])), f"Unique identifier should be {chrom_type} + 0 for paternal, and {chrom_type} + 23 for maternal (entry: {entry})"
            assert entry["parent"] == -1, f"Parent of root node should be -1 (entry: {entry})"
            assert entry["epoch_created"] == 0, f"Epoch created should be 0 (entry: {entry})"
            assert isinstance(entry["paternal"], bool), f"Paternal attribute should be a boolean value (entry: {entry})"
            assert entry["SNVs"] == [], f"SNVs should be initialized as an empty list (entry: {entry})"
            assert entry["dead"] == False, f"Dead attribute should be False initially (entry: {entry})"


test_initialize_simulated_chromosomes()


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


def test_simulate_snvs():
    chromosomes = initialize_simulated_chromosomes()
    rate = 0.05
    epoch = 1

    initial_snv_count = 0
    for chrom_type, chrom_data in chromosomes.items():
        for entry in chrom_data:
            initial_snv_count += len(entry["SNVs"])

    final_snv_count, updated_chromosomes = simulate_snvs(chromosomes, LENGTHS, rate, epoch, initial_snv_count)
    assert final_snv_count >= initial_snv_count, "SNV count should increase after simulation"

    for chrom_type, chrom_data in updated_chromosomes.items():
        for entry in chrom_data:
            SNVS = entry["SNVs"]
            assert type(SNVS) == list
            for snv in SNVS:
                assert type(snv) == dict
                assert 'unique_identifier' in snv
                assert type(snv['unique_identifier']) == str
                assert 'epoch_created' in snv
                assert type(snv['epoch_created']) == int
                assert snv['epoch_created'] == 1


test_simulate_snvs()

from typing import Dict, List

def simulate_gd(
    simulated_chromosomes: Dict[str, List[Dict]], 
    epoch: int, 
    chrom_count: int
    ) -> (int, Dict[str, List[Dict]]):
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


def create_new_chromosome(
    old_chromosome: Dict, 
    chrom_count: int, 
    epoch: int
    ) -> Dict:
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


def test_simulate_gd():
    chromosomes = initialize_simulated_chromosomes()
    epoch = 1
    initial_chrom_count = 46

    final_chrom_count, updated_chromosomes = simulate_gd(chromosomes, epoch, initial_chrom_count)
    assert final_chrom_count == initial_chrom_count * 2, "Chromosome count should double after genome doubling"

    for chrom_type, chrom_data in updated_chromosomes.items():
        for entry in chrom_data:
            if entry["epoch_created"] == epoch:
                assert entry["parent"] != -1, "Parent should not be -1 for newly created chromosomes"

    # Count the number of chromosomes in updated_chromosomes
    actual_chrom_count = sum(len(chrom_list) for chrom_list in updated_chromosomes.values())

    # Check if every chrom_type has an even number of chromosomes
    for chrom_type, chrom_list in updated_chromosomes.items():
        assert len(chrom_list) % 2 == 0, f"Chrom_type {chrom_type} should have an even number of chromosomes"

    # Check if actual chromosome count in updated_chromosomes is even
    assert actual_chrom_count % 2 == 0, "Actual chromosome count in updated_chromosomes should be an even number"

    # Check if actual chromosome count in updated_chromosomes matches the final_chrom_count
    assert actual_chrom_count == final_chrom_count, "Actual chromosome count in updated_chromosomes should match the final_chrom_count"

    # Check if half of each chromosome type have epoch = 2
    for chrom_type, chrom_list in updated_chromosomes.items():
        total_chromosomes_type = len(chrom_list)
        chromosomes_epoch_2_type = sum(chrom['epoch_created'] == 1 for chrom in chrom_list)

        assert chromosomes_epoch_2_type == total_chromosomes_type / 2, f"Half of the chromosomes of type {chrom_type} should have epoch = 2"


test_simulate_gd()

def simulate_anueploidy_agnostic(
    simulated_chromosomes: Dict[str, List[Dict]], 
    epoch: int, 
    chrom_count: int
) -> (int, Dict[str, List[Dict]]):
    """
    Simulate Anueploidy on a given set of chromosomes, agnostic of chromosome details.

    Randomly select chromosomes to lose and gain, updating their status accordingly.

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
        chrom_count = simulate_agnostic_chromosome_loss(chromosomes, chrom_count)
        chrom_count = simulate_agnostic_chromosome_gain(chromosomes, epoch, chrom_count)

    return chrom_count, simulated_chromosomes


def simulate_agnostic_chromosome_loss(chromosomes: List[Dict], chrom_count: int) -> int:
    """
    Simulate the loss of a random number of chromosomes.

    Args:
        chromosomes (List[Dict]): A list of chromosomes represented as dictionaries.

        chrom_count (int): The current chromosome count.

    Returns:
        chrom_count (int): The updated chromosome count after the loss.
    """
    for which in random.sample(range(len(chromosomes)), random.randint(0, len(chromosomes))):
        chromosomes[which]["dead"] = True
        chrom_count -= 1

    return chrom_count


def simulate_agnostic_chromosome_gain(chromosomes: List[Dict], epoch: int, chrom_count: int) -> int:
    """
    Simulate the gain of a random number of chromosomes.

    Args:
        chromosomes (List[Dict]): A list of chromosomes represented as dictionaries.

        epoch (int): The current epoch.

        chrom_count (int): The current chromosome count.

    Returns:
        chrom_count (int): The updated chromosome count after the gain.
    """
    for which in random.choices(range(len(chromosomes)), k=random.randint(0, len(chromosomes))):
        if not chromosomes[which]["dead"]:
            new_chrom = create_new_chromosome(chromosomes[which], chrom_count, epoch)
            chromosomes.append(new_chrom)
            chrom_count += 1

    return chrom_count



def test_simulate_anueploidy_agnostic():
    chromosomes = initialize_simulated_chromosomes()
    epoch = 1
    chrom_count = 46

    updated_chrom_count, updated_chromosomes = simulate_anueploidy_agnostic(chromosomes, epoch, chrom_count)
    assert 0 <= updated_chrom_count <= 2 * chrom_count, "Chromosome count should be within the range [0, 2 * chrom_count] after 1 round of simulation"


    # Check all newly created chromosomes have epoch == 1
    for chrom_list in updated_chromosomes.values():
        for chrom in chrom_list:
            if chrom["epoch_created"] > epoch - 1:  # for newly created chromosomes
                assert chrom["epoch_created"] == 1, "Newly created chromosomes should have epoch_created = 1"


test_simulate_anueploidy_agnostic()


def simulate_anueploidy(
    simulated_chromosomes: Dict[str, List[Dict]], 
    epoch: int, 
    chrom_count: int, 
    p_up: float, 
    p_down: float
    ) -> (int, Dict[str, List[Dict]]):
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


def simulate_chromosome_changes(
    chromosomes: List[Dict], 
    epoch: int, 
    chrom_count: int, 
    p_up: float, 
    p_down: float
    ) -> List[Dict]:
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


def test_simulate_anueploidy():
    chromosomes = initialize_simulated_chromosomes()
    epoch = 1
    chrom_count = 46
    p_up = 0.1
    p_down = 0.1

    updated_chrom_count, updated_chromosomes = simulate_anueploidy(chromosomes, epoch, chrom_count, p_up, p_down)
    assert updated_chrom_count >= chrom_count, "Chromosome count should increase or stay the same after simulation"


    # Check all newly created chromosomes have epoch == 1
    for chrom_list in updated_chromosomes.values():
        for chrom in chrom_list:
            if chrom["epoch_created"] > epoch - 1:  # for newly created chromosomes
                assert chrom["epoch_created"] == 1, "Newly created chromosomes should have epoch_created = 1"

test_simulate_anueploidy()


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


def simulate_cancer_genome(p_up: float, p_down: float, pre: int, mid: int, post: int, rate: float, agnostic: bool=False) -> List:
    """
    Simulate a genome based on given parameters.

    :param p_up: Probability of upward mutation.
    :param p_down: Probability of downward mutation.
    :param pre: Pre-mutation period.
    :param mid: Mid-mutation period.
    :param post: Post-mutation period.
    :param rate: Mutation rate.
    :param agnostic: Indicates if the simulation is agnostic or not.
    :return: Simulated chromosomes.
    """
    # Count of unique SNVs
    snv_count = 0

    # Initialize chromosomes
    simulated_chromosomes = initialize_simulated_chromosomes()

    # Count of unique chromosomes
    chrom_count = 46  # in a human cell

    # Total epochs count
    ev_sequence = get_ev_string(pre, mid, post)

    if len(ev_sequence) == 0:
        return simulated_chromosomes

    for epoch, epoch_type in enumerate(ev_sequence):
        # Simulate SNVs
        # note that SNVs are created in this epoch, whilst we "time" the copy number changes to be in the next epoch when they are created.
        snv_count, simulated_chromosomes = simulate_snvs(simulated_chromosomes, LENGTHS, rate, epoch, snv_count)
        check_all_chrs_are_unique(simulated_chromosomes)

        if (mid != -1 and epoch == pre) or (post != -1 and epoch == pre + 1 + mid):
            assert epoch_type == "G"
            chrom_count, simulated_chromosomes = simulate_gd(simulated_chromosomes, epoch+1, chrom_count)
            check_all_chrs_are_unique(simulated_chromosomes)
        else:
            assert epoch_type == "A"
            if agnostic:
                chrom_count, simulated_chromosomes = simulate_anueploidy_agnostic(simulated_chromosomes, epoch+1, chrom_count)
                check_all_chrs_are_unique(simulated_chromosomes)
            else:
                chrom_count, simulated_chromosomes = simulate_anueploidy(simulated_chromosomes, epoch+1, chrom_count, p_up, p_down)
                check_all_chrs_are_unique(simulated_chromosomes)

    assert pre + mid + post + 2 == len(ev_sequence)
    check_simulated_chromosomes(simulated_chromosomes, pre, mid, post, ev_sequence)

    return simulated_chromosomes



def test_simulate_cancer_genome():
    p_up, p_down = 0.1, 0.1
    pre, mid, post = 2, 1, 3
    rate = 0.01
    agnostic = False

    if mid == -1:
        assert post == -1, "If mid is -1, post should also be -1"

    simulated_chromosomes = simulate_cancer_genome(
        p_up, p_down, pre, mid, post, rate, agnostic
    )

    for chrom_type in simulated_chromosomes:
        assert len(simulated_chromosomes) != 0, "There should be chromosomes in the simulation"

    if post == 0 or (mid == 0 and post == -1):
        assert simulated_chromosomes[-1] == "G", "If a genome doubling round just occurred, the last element of the simulated_chromosomes should be 'G'"

    assert check_all_chrs_are_unique(simulated_chromosomes), "All chromosomes in the simulation should be unique"
    assert check_expected_keys_in_simulated_chromosomes_present(simulated_chromosomes)


test_simulate_cancer_genome()


##### STEP 1b; from the simulated genome create a tree
#####
#####
#####
#####
#####


def insert_node_into_leaf(tree: Dict, node: Dict) -> None:
    """
    Inserts a node into a leaf of a given tree.

    :param tree: The tree where the node is to be inserted.
    :param node: The node to be inserted into the tree.
    """

    # Assert that the given tree is a leaf node
    assert tree["child"] is None, f"Expected tree['child'] to be None, but got {tree['child']}"
    assert tree["complement"] is None, f"Expected tree['complement'] to be None, but got {tree['complement']}"

    # Create a deep copy of the tree as the complement
    tree["complement"] = copy.deepcopy(tree)
    tree["complement"]["epoch_created"] = node["epoch_created"]

    # Insert the node as a child
    tree["child"] = copy.deepcopy(node)
    tree["child"]["child"] = None
    tree["child"]["complement"] = None

def test_insert_node_into_leaf():
    # Case 1: Insert a node under the child
    tree = {
        "unique_identifier": 1,
        "child": None,
        "complement": None
    }
    node = {
        "unique_identifier": 2,
        "parent": 1,
        "epoch_created": 100
    }
    expected_tree = {
        "unique_identifier": 1,
        "child": {
            "unique_identifier": 2,
            "parent": 1,
            "epoch_created": 100,
            "child": None,
            "complement": None
        },
        "complement": {
            "unique_identifier": 1,
            "child": None,
            "complement": None,
            "epoch_created": 100
        }
    }
    insert_node_into_leaf(tree, node)
    assert tree == expected_tree, f"Expected {expected_tree}, but got {tree}"


def insert_node_under_complement(tree: Dict, node: Dict) -> None:
    """
    Inserts a node under the complement of a given tree.

    :param tree: The tree where the node is to be inserted.
    :param node: The node to be inserted under the complement of the tree.
    """

    # nodes are inserted into the tree in their order of their unique identifier
    # which maps to the order of chromosomes created.
    # there may be many children of "complement" which is why it gets dragged down the tree
    # to find out when it bifurcated into each new chromosome it is the parent of

    # complement cannot be None, because the "child" node is always created first
    # and complement is always created with it as a place holder
    assert tree["complement"] is not None, f"Expected tree['complement'] to not be None, but it was."

    if node["parent"] == -1:
        # these can't get passed further on because they are the original chromosomes
        expected_keys = {"unique_identifier", "epoch_created", "parent", "SNVs", "paternal", "dead"}
        assert set(node.keys()) == expected_keys, f"Unexpected keys: {node.keys()} - Expected keys: {expected_keys}"

        assert node["unique_identifier"] < 46, f"Expected node['unique_identifier'] to be less than 46, but got {node['unique_identifier']}"

        tree["complement"] = copy.deepcopy(node)
        tree["complement"]["child"] = None
        tree["complement"]["complement"] = None
    else:
        # we don't know how far down the complement line that we have to insert the node: so just keep going until it is done:

        # insert_node_under_complement and insert_node_under_leaf do not return anything,
        # they modify the tree structure given to them,
        # why do we return something here from insert_node_into_tree_structure? CHECKHERE

        tree["complement"] = insert_node_into_tree_structure(tree["complement"], node)


def test_insert_node_under_complement():
    # Case 1: Insert a node under the complement
    tree = {
        "unique_identifier": 1,
        "child": {
            "unique_identifier": 2,
            "parent": 1,
            "child": None,
            "complement": None
        },
        "complement": {
            "unique_identifier": 1,
            "child": None,
            "complement": None
        }
    }
    node = {
        "unique_identifier": 3,
        "parent": 1,
        "epoch_created": 1
    }
    expected_tree = {
        "unique_identifier": 1,
        "child": {
            "unique_identifier": 2,
            "parent": 1,
            "child": None,
            "complement": None
        },
        "complement": {
            "unique_identifier": 1,
            "child": {
                "unique_identifier": 3,
                "parent": 1,
                "epoch_created": 1,
                "child": None,
                "complement": None
            },
            "complement": {
                "unique_identifier": 1,
                "child": None,
                "complement": None,
                "epoch_created": 1
            }
        }
    }
    insert_node_under_complement(tree, node)
    assert tree == expected_tree, f"Expected:\n{json.dumps(expected_tree, indent=4)}\nBut got:\n{json.dumps(tree, indent=4)}"


def insert_node_into_tree_structure(tree: Optional[Dict], node: Optional[Dict]) -> Optional[Dict]:
    """
    Inserts a node into a tree structure according to the unique identifier and parent-child relationships.

    :param tree: The tree where the node is to be inserted.
    :param node: The node to be inserted into the tree.
    :return: The modified tree.
    """

    if node is None or tree is None:
        return tree

    if node["unique_identifier"] is not None and tree["unique_identifier"] is not None and node["unique_identifier"] != tree["unique_identifier"]:
        if node["parent"] == tree["unique_identifier"]:
            if tree["child"] is None:
                insert_node_into_leaf(tree, node)
            else:
                insert_node_under_complement(tree, node)
        else:
            if tree["child"] is not None:
                tree["child"] = insert_node_into_tree_structure(tree["child"], node)

            if tree["complement"] is not None:
                tree["complement"] = insert_node_into_tree_structure(tree["complement"], node)

    return tree



def test_insert_node_into_tree_structure():
    # Test case 1
    tree1 = {
        "unique_identifier": 1,
        "child": None,
        "complement": None,
        "epoch_created": 1
    }
    node1 = {
        "unique_identifier": 2,
        "parent": 1,
        "epoch_created": 2
    }
    expected_tree1 = {
        "unique_identifier": 1,
        "child": {
            "unique_identifier": 2,
            "parent": 1,
            "epoch_created": 2,
            "child": None,
            "complement": None
        },
        "complement": {
            "unique_identifier": 1,
            "child": None,
            "complement": None,
            "epoch_created": 2
        },
        "epoch_created": 1
    }
    result_tree1 = insert_node_into_tree_structure(tree1, node1)
    assert result_tree1 == expected_tree1, f"Test case 1 failed: Expected:\n{json.dumps(expected_tree1, indent=4)}\nBut got:\n{json.dumps(result_tree1, indent=4)}"

    # Test case 2
    tree2 = {
        "unique_identifier": 1,
        "child": {
            "unique_identifier": 2,
            "parent": 1,
            "child": None,
            "complement": None,
            "epoch_created": 2
        },
        "complement": None,
        "epoch_created": 1
    }
    node2 = {
        "unique_identifier": 3,
        "parent": 2,
        "epoch_created": 3
    }
    expected_tree2 = {
        "unique_identifier": 1,
        "child": {
            "unique_identifier": 2,
            "parent": 1,
            "child": {
                "unique_identifier": 3,
                "parent": 2,
                "epoch_created": 3,
                "child": None,
                "complement": None
            },
            "complement": {
                "unique_identifier": 2,
                "parent": 1,
                "child": None,
                "complement": None,
                "epoch_created": 3
            },
            "epoch_created": 2
        },
        "complement": None,
        "epoch_created": 1
    }
    result_tree2 = insert_node_into_tree_structure(tree2, node2)
    assert result_tree2 == expected_tree2, f"Test case 2 failed: Expected:\n{json.dumps(expected_tree2, indent=4)}\nBut got:\n{json.dumps(result_tree2, indent=4)}"

    # Test case 3
    tree3 = None
    node3 = {
        "unique_identifier": 1,
        "parent": None,
        "epoch_created": 2
    }
    expected_tree3 = None
    result_tree3 = insert_node_into_tree_structure(tree3, node3)
    assert result_tree3 == expected_tree3, f"Test case 3 failed: Expected:\n{json.dumps(expected_tree3, indent=4)}\nBut got:\n{json.dumps(result_tree3, indent=4)}"

    # Test case 4
    tree4 = {
        "unique_identifier": 1,
        "child": None,
        "complement": None
    }
    node4 = None
    expected_tree4 = {
        "unique_identifier": 1,
        "child": None,
        "complement": None
    }
    result_tree4 = insert_node_into_tree_structure(tree4, node4)
    assert result_tree4 == expected_tree4, f"Test case 4 failed: Expected:\n{json.dumps(expected_tree4, indent=4)}\nBut got:\n{json.dumps(result_tree4, indent=4)}"


test_insert_node_into_leaf()
test_insert_node_under_complement()
test_insert_node_into_tree_structure()


def add_copynumber_to_tree_structure(tree: Optional[Dict]) -> Optional[Dict]:
    """
    Add the "copy_number" attribute to each node in the given tree structure.

    :param tree: The tree to which the "copy_number" attribute is to be added.
    :return: The modified tree with added "copy_number" attribute.
    """

    if "copy_number" in tree:
        return tree

    if tree["child"] is None and tree["complement"] is None:
        tree["copy_number"] = 0 if tree.get("dead") else 1
    else:
        if tree["child"] is not None:
            tree["child"] = add_copynumber_to_tree_structure(tree["child"])
        if tree["complement"] is not None:
            tree["complement"] = add_copynumber_to_tree_structure(tree["complement"])

        child_copy_number = tree["child"]["copy_number"] if tree["child"] else 0
        complement_copy_number = tree["complement"]["copy_number"] if tree["complement"] else 0
        tree["copy_number"] = child_copy_number + complement_copy_number

    return tree



def test_add_copynumber_to_tree_structure():
    # Test case 1: Simple tree with no children or complements
    tree = {
        "child": None,
        "complement": None,
        "dead": False,
    }
    expected_tree = {
        "child": None,
        "complement": None,
        "dead": False,
        "copy_number": 1,
    }
    result_tree = add_copynumber_to_tree_structure(tree)
    assert result_tree == expected_tree, f"Expected {expected_tree}, but got {result_tree}"

    # Test case 2: Tree with one level of children
    tree = {
        "child": {
            "child": None,
            "complement": None,
            "dead": False,
        },
        "complement": {
            "child": None,
            "complement": None,
            "dead": True,
        },
    }
    expected_tree = {
        "child": {
            "child": None,
            "complement": None,
            "dead": False,
            "copy_number": 1,
        },
        "complement": {
            "child": None,
            "complement": None,
            "dead": True,
            "copy_number": 0,
        },
        "copy_number": 1,
    }
    result_tree = add_copynumber_to_tree_structure(tree)
    assert result_tree == expected_tree, f"Expected {expected_tree}, but got {result_tree}"

    # Test case 3: Tree with multiple levels
    tree = {
        "child": {
            "child": {
                "child": None,
                "complement": None,
                "dead": False,
            },
            "complement": {
                "child": None,
                "complement": None,
                "dead": True,
            },
        },
        "complement": {
            "child": None,
            "complement": None,
            "dead": False,
        },
    }
    expected_tree = {
        "child": {
            "child": {
                "child": None,
                "complement": None,
                "dead": False,
                "copy_number": 1,
            },
            "complement": {
                "child": None,
                "complement": None,
                "dead": True,
                "copy_number": 0,
            },
            "copy_number": 1,
        },
        "complement": {
            "child": None,
            "complement": None,
            "dead": False,
            "copy_number": 1,
        },
        "copy_number": 2,
    }
    result_tree = add_copynumber_to_tree_structure(tree)
    assert result_tree == expected_tree, f"Expected {expected_tree}, but got {result_tree}"

test_add_copynumber_to_tree_structure()


def add_SNV_multiplicity_to_tree_structure(tree: Optional[Dict]) -> Optional[Dict]:
    """
    Adds "SNV_multiplicity" attribute to each node in the given tree structure.

    :param tree: The tree to which the "SNV_multiplicity" attribute is to be added.
    :return: The modified tree with added "SNV_multiplicity" attribute.
    """

    if tree["copy_number"] == 0:
        tree["SNV_multiplicity"] = None
        return tree

    if tree["child"] is None:
        assert tree["complement"] is None
        count = sum(1 for SNV in tree["SNVs"] if SNV["epoch_created"] >= tree["epoch_created"])
        tree["SNV_multiplicity"] = count
        return tree

    assert tree["child"]["epoch_created"] == tree["complement"]["epoch_created"]

    tree["child"] = add_SNV_multiplicity_to_tree_structure(tree["child"])
    tree["complement"] = add_SNV_multiplicity_to_tree_structure(tree["complement"])

    count = sum(1 for SNV in tree["SNVs"] if tree["epoch_created"] <= SNV["epoch_created"] < tree["child"]["epoch_created"])
    tree["SNV_multiplicity"] = count

    return tree


def test_add_SNV_multiplicity_to_tree_structure():
    # Test case 1: simple tree with no children or complements
    tree = {
        "copy_number": 1,
        "child": None,
        "complement": None,
        "epoch_created": 0,
        "SNVs": [{"epoch_created": 0}, {"epoch_created": 1}]
    }
    expected_tree = tree.copy()
    expected_tree["SNV_multiplicity"] = 2
    result_tree = add_SNV_multiplicity_to_tree_structure(tree)
    assert result_tree == expected_tree, f"Expected {expected_tree}, but got {result_tree}"

    # Test case 2: tree with one level of children
    tree = {
        "copy_number": 2,
        "epoch_created": 0,
        "SNVs": [{"epoch_created": 0}],
        "child": {
            "copy_number": 1,
            "epoch_created": 1,
            "SNVs": [{"epoch_created": 1}],
            "child": None,
            "complement": None
        },
        "complement": {
            "copy_number": 1,
            "epoch_created": 1,
            "SNVs": [{"epoch_created": 2}],
            "child": None,
            "complement": None
        }
    }
    expected_tree = {
        "copy_number": 2,
        "epoch_created": 0,
        "SNV_multiplicity": 1,
        "SNVs": [{"epoch_created": 0}],
        "child": {
            "copy_number": 1,
            "epoch_created": 1,
            "SNV_multiplicity": 1,
            "SNVs": [{"epoch_created": 1}],
            "child": None,
            "complement": None
        },
        "complement": {
            "copy_number": 1,
            "epoch_created": 1,
            "SNV_multiplicity": 1,
            "SNVs": [{"epoch_created": 2}],
            "child": None,
            "complement": None
        }
    }
    result_tree = add_SNV_multiplicity_to_tree_structure(tree)
    assert result_tree == expected_tree, f"Expected {expected_tree}, but got {result_tree}"

    # Additional test cases can be added for more complex trees

test_add_SNV_multiplicity_to_tree_structure()


def remove_SNVs_from_tree_structure(tree: Optional[Dict]) -> Optional[Dict]:
    """
    Recursively removes the "SNVs" attribute from each node in the given tree structure.

    :param tree: The tree from which the "SNVs" attribute is to be removed.
    :return: The modified tree without the "SNVs" attribute.
    """

    tree.pop("SNVs", None)

    if tree["child"] is not None:
        assert tree["complement"] is not None
        tree["child"] = remove_SNVs_from_tree_structure(tree["child"])
        tree["complement"] = remove_SNVs_from_tree_structure(tree["complement"])

    return tree



def test_remove_SNVs_from_tree_structure():
    # Test case 1
    tree = {
        "SNVs": ["dummy"],
        "child": None,
        "complement": None
    }
    expected_tree = {
        "child": None,
        "complement": None
    }
    result_tree = remove_SNVs_from_tree_structure(tree)
    assert result_tree == expected_tree, f"Expected {expected_tree}, but got {result_tree}"


def remove_dead_leaf_nodes_from_tree_structure(tree: Optional[Dict]) -> Optional[Dict]:
    """
    Recursively removes dead leaf nodes from the given tree structure.

    :param tree: The tree from which dead leaf nodes are to be removed.
    :return: The modified tree with dead leaf nodes removed.
    """

    if tree is None:
        return None

    tree['child'] = remove_dead_leaf_nodes_from_tree_structure(tree['child'])
    tree['complement'] = remove_dead_leaf_nodes_from_tree_structure(tree['complement'])

    if tree['child'] is None and tree['complement'] is None and tree.get('dead', False) and tree['parent'] != -1:
        return None

    return tree


def test_remove_dead_leaf_nodes_from_tree_structure():
    # Test case 1
    tree = {
        "dead": True,
        "parent": 1,
        "child": None,
        "complement": None
    }
    expected_tree = None
    result_tree = remove_dead_leaf_nodes_from_tree_structure(tree)
    assert result_tree == expected_tree, f"Expected {expected_tree}, but got {result_tree}"



def are_all_descendants_zero_in_tree_structure(node: Optional[Dict]) -> bool:
    """
    Checks if the "copy_number" attribute of the given node and all its descendants in the tree structure is zero.

    :param node: The node to start the check from.
    :return: True if all "copy_number" attributes are zero, False otherwise.
    """

    if node is None:
        return True

    if node["copy_number"] != 0:
        return False

    return are_all_descendants_zero_in_tree_structure(node["child"]) and are_all_descendants_zero_in_tree_structure(node["complement"])


def test_are_all_descendants_zero_in_tree_structure():
    # Test case 1
    node = {
        "copy_number": 0,
        "child": None,
        "complement": None
    }
    expected_result = True
    result = are_all_descendants_zero_in_tree_structure(node)
    assert result == expected_result, f"Expected {expected_result}, but got {result}"

    # Add more test cases as needed

test_remove_SNVs_from_tree_structure()
test_remove_dead_leaf_nodes_from_tree_structure()
test_are_all_descendants_zero_in_tree_structure()


def check_and_remove_redundant_node_in_tree_structure(node: Optional[Dict], main_child: Optional[Dict], other_child: Optional[Dict]) -> Optional[Dict]:
    """
    Checks if a node in the tree structure is redundant based on its children's "copy_number" attributes,
    and if it is, it updates the "epoch_created" of the main child and removes the node.

    :param node: The node to check.
    :param main_child: The main child node of the node.
    :param other_child: The other child node of the node.
    :return: The updated node (if it was redundant, it gets replaced with the main child node).
    """

    if main_child is not None and main_child["copy_number"] == node["copy_number"]:
        if not are_all_descendants_zero_in_tree_structure(other_child):
            raise ValueError("Invalid tree: child and complement copy_number do not add up to parent's copy_number")

        # Update epoch_created of the main_child to that of its parent
        main_child["epoch_created"] = node["epoch_created"]
        return main_child

    return node


def remove_redundant_parents_from_tree_structure(node: Optional[Dict]) -> Optional[Dict]:
    """
    Recursively removes redundant parents from a tree structure. A node is considered redundant if its "copy_number"
    equals its child's "copy_number" and all the descendants of the other child node have a "copy_number" of zero.

    :param node: The root of the tree.
    :return: The root of the tree after redundant nodes have been removed.
    """

    if node is None:
        return None
        
    # Recursively process child and complement nodes
    node["child"] = remove_redundant_parents_from_tree_structure(node["child"])
    node["complement"] = remove_redundant_parents_from_tree_structure(node["complement"])
    
    # Check for the conditions to remove the node, only if the node has a parent
    if node["parent"] is not None:
        node = check_and_remove_redundant_node_in_tree_structure(node, node["child"], node["complement"])
        node = check_and_remove_redundant_node_in_tree_structure(node, node["complement"], node["child"])
        
    return node


def test_remove_redundant_parents_from_tree_structure():
    # Test case 1
    node = {
        "parent": 1,
        "copy_number": 1,
        "child": {
            "parent": 1,
            "copy_number": 1,
            "child": None,
            "complement": None,
            "epoch_created": 2
        },
        "complement": {
            "parent": 1,
            "copy_number": 0,
            "child": None,
            "complement": None,
            "epoch_created": 2
        },
        "epoch_created": 1
    }
    expected_node = {
        "parent": 1,
        "copy_number": 1,
        "child": None,
        "complement": None,
        "epoch_created": 1
    }
    result_node = remove_redundant_parents_from_tree_structure(node)
    assert result_node == expected_node, f"Expected {expected_node}, but got {result_node}"

    # Add more test cases as needed

test_remove_redundant_parents_from_tree_structure()


def assign_epoch_killed_to_tree_structure(tree, max_epochs):
    """
    The function calculates the epoch_killed for each node in the given tree.

    :param tree: A dictionary representing a tree node with keys:
                 'child', 'complement', 'epoch_created', and 'epoch_killed'.
    :return: The modified tree with 'epoch_killed' values added.
    """

    if tree is None:
        return None

    child = tree.get('child')
    complement = tree.get('complement')

    if child is not None or complement is not None:
        child_epoch_created = child.get('epoch_created')
        complement_epoch_created = complement.get('epoch_created')

        if child_epoch_created != complement_epoch_created:
            raise ValueError("Epoch created values of child and complement do not match.")

        tree['child'] = assign_epoch_killed_to_tree_structure(child, max_epochs)
        tree['complement'] = assign_epoch_killed_to_tree_structure(complement, max_epochs)

        tree['epoch_killed'] = child['epoch_created']
    else:
        tree['epoch_killed'] = max_epochs

    return tree

def test_assign_epoch_killed_to_tree_structure():
    # Test case 1
    tree = {
        "epoch_created": 1,
        "child": {
            "epoch_created": 2,
            "child": None,
            "complement": None
        },
        "complement": {
            "epoch_created": 2,
            "child": None,
            "complement": None
        }
    }
    max_epochs = 3
    expected_tree = {
        "epoch_created": 1,
        "epoch_killed": 2,
        "child": {
            "epoch_created": 2,
            "epoch_killed": 3,
            "child": None,
            "complement": None
        },
        "complement": {
            "epoch_created": 2,
            "epoch_killed": 3,
            "child": None,
            "complement": None
        }
    }
    result_tree = assign_epoch_killed_to_tree_structure(tree, max_epochs)
    assert result_tree == expected_tree, f"Expected {expected_tree}, but got {result_tree}"

    # Add more test cases as needed

test_assign_epoch_killed_to_tree_structure()



def create_truth_trees_from_simulation(simulated_chromosomes: Dict, max_epochs: int) -> Dict:
    """
    Given the simulated chromosomes, the function creates the truth trees.
    It creates a separate tree for each type of chromosome.
    Each tree is updated and simplified by adding metadata to the nodes,
    adding copy numbers, SNV multiplicity, removing SNVs and dead nodes, 
    removing redundant nodes, and assigning 'epoch_killed' to nodes.

    :param simulated_chromosomes: A dictionary of simulated chromosomes.
    :param max_epochs: The maximum number of epochs for simulation.
    :return: A dictionary of truth trees for each chromosome type.
    """
    # There is a tree for every chromosome
    trees = {}

    for chrom_type in simulated_chromosomes:
        # First sort the nodes by the order they need to be inserted in:
        sorted_list = sorted([(x["unique_identifier"],x) for x in simulated_chromosomes[chrom_type]])

        # Create the root node of the tree for this chrom_type
        tree = {
            'unique_identifier': -1,
            'parent': None,
            'epoch_created': 0,
            'paternal': None,
            'child': None,
            'complement': None,
            'SNVs': []
        }

        # Insert all nodes and add metadata to tree and simplify:
        for new_node in sorted_list:
            trees[chrom_type] = insert_node_into_tree_structure(tree, new_node[1])  # UNSURE EDIT

        logging.debug("##### chrom_type: %s", str(chrom_type))
        logging.debug("With nodes inserted:")
        logging.debug("%s", trees[chrom_type])

        trees[chrom_type] = add_copynumber_to_tree_structure(trees[chrom_type])
        logging.debug("With copynumber annotated:")
        logging.debug("%s", trees[chrom_type])

        trees[chrom_type] = add_SNV_multiplicity_to_tree_structure(trees[chrom_type])
        logging.debug("With SNV multiplicity:")
        logging.debug("%s", trees[chrom_type])

        trees[chrom_type] = remove_SNVs_from_tree_structure(trees[chrom_type])
        logging.debug("With SNVs removed from the tree:")
        logging.debug("%s", trees[chrom_type])

        trees[chrom_type] = remove_dead_leaf_nodes_from_tree_structure(trees[chrom_type])
        logging.debug("With dead nodes removed:")
        logging.debug("%s", trees[chrom_type])

        trees[chrom_type] = remove_redundant_parents_from_tree_structure(trees[chrom_type])
        logging.debug("With redundant nodes removed:")
        logging.debug("%s", trees[chrom_type])

        trees[chrom_type] = assign_epoch_killed_to_tree_structure(trees[chrom_type], max_epochs)
        logging.debug("With epoch killed:")
        logging.debug("%s", trees[chrom_type])

    return trees


def get_child_and_complement_trees_from_truth_tree(truth_tree: Dict) -> Tuple:
    """
    Get child and complement trees from a given truth tree.

    :param truth_tree: The input truth tree.
    :return: The child and complement trees.
    """
    child_tree = convert_truth_tree_to_CN_tree(truth_tree["child"]) if truth_tree["child"] else None
    complement_tree = convert_truth_tree_to_CN_tree(truth_tree["complement"]) if truth_tree["complement"] else None
    return child_tree, complement_tree


def convert_truth_tree_to_CN_tree(truth_tree: Dict) -> List:
    """
    Convert truth tree to copy number (CN) tree.

    :param truth_tree: The input truth tree.
    :return: The resultant CN tree.
    """
    child_tree, complement_tree = get_child_and_complement_trees_from_truth_tree(truth_tree)

    if child_tree and complement_tree:
        CN_tree = [truth_tree["copy_number"], child_tree, complement_tree]
    elif child_tree:
        CN_tree = [truth_tree["copy_number"], child_tree]
    elif complement_tree:
        CN_tree = [truth_tree["copy_number"], complement_tree]
    else:
        CN_tree = [truth_tree["copy_number"]]

    return CN_tree

def test_convert_truth_tree_to_CN_tree():
    # Case 1: Test with a simple truth_tree
    truth_tree = {
        "copy_number": 1,
        "child": {
            "copy_number": 2,
            "child": None,
            "complement": None
        },
        "complement": {
            "copy_number": 3,
            "child": None,
            "complement": None
        }
    }
    expected_CN_tree = [1, [2], [3]]
    result_CN_tree = convert_truth_tree_to_CN_tree(truth_tree)
    assert result_CN_tree == expected_CN_tree, f"Expected {expected_CN_tree}, but got {result_CN_tree}"

test_convert_truth_tree_to_CN_tree()


def convert_truth_trees_to_CN_trees(truth_trees: Dict[str, Dict]) -> Dict[str, List]:
    """
    Convert truth trees to copy number (CN) trees and make the trees left-heavy.

    :param truth_trees: The input truth trees.
    :return: The resultant CN trees.
    """
    for chrom_type in truth_trees:
        truth_trees[chrom_type] = convert_truth_tree_to_CN_tree(truth_trees[chrom_type])
        truth_trees[chrom_type] = make_left_heavy(truth_trees[chrom_type])
    return truth_trees


def test_convert_truth_trees_to_CN_trees():
    # Case 1: Test with simple truth_trees
    truth_trees = {
        "type1": {
            "copy_number": 1,
            "child": {
                "copy_number": 2,
                "child": None,
                "complement": None
            },
            "complement": {
                "copy_number": 3,
                "child": None,
                "complement": None
            }
        }
    }
    expected_CN_trees = {
        "type1": [1, [3], [2]]
    }
    result_CN_trees = convert_truth_trees_to_CN_trees(truth_trees)
    assert result_CN_trees == expected_CN_trees, f"Expected {expected_CN_trees}, but got {result_CN_trees}"

test_convert_truth_trees_to_CN_trees()


def count_copy_numbers(simulated_chromosomes: Dict[str, List[Dict]]) -> Dict[str, List[int]]:
    """
    Count the number of each parental specific copy number found in the genome.

    :param simulated_chromosomes: The input simulated chromosomes.
    :return: A dictionary mapping chromosome type to a list of counts of living chromosomes for each parent.
    """
    observed_copy_numbers = {}
    for chrom_type in simulated_chromosomes:
        observed_copy_numbers[chrom_type] = [
            len([x for x in simulated_chromosomes[chrom_type] if paternal == x["paternal"] and not x["dead"]])
            for paternal in [True, False]
        ]
    return observed_copy_numbers

def test_count_copy_numbers():
    simulated_chromosomes = {
        "type1": [
            {"paternal": True, "dead": False},
            {"paternal": True, "dead": True},
            {"paternal": False, "dead": False},
            {"paternal": False, "dead": False},
        ],
        "type2": [
            {"paternal": True, "dead": False},
            {"paternal": True, "dead": False},
            {"paternal": False, "dead": True},
            {"paternal": False, "dead": True},
        ],
    }

    expected_counts = {
        "type1": [1, 2],  # 1 living paternal and 2 living non-paternal for type1
        "type2": [2, 0],  # 2 living paternal and 0 living non-paternal for type2
    }

    observed_counts = count_copy_numbers(simulated_chromosomes)

    assert observed_counts == expected_counts, f"Expected {expected_counts}, but got {observed_counts}"

test_count_copy_numbers()


def count_CN_multiplicities(observed_CNs):
    multiplicities = {}

    for chrom_type in observed_CNs:
        for CN in observed_CNs[chrom_type]:
            if CN not in multiplicities:
                multiplicities[CN] = 1

            else:
                multiplicities[CN] += 1

    return(multiplicities)


def count_copy_number_multiplicities(observed_copy_numbers: Dict[str, List[int]]) -> Dict[int, int]:
    """
    Count the multiplicities of each observed copy number in the genome.

    :param observed_copy_numbers: The input observed copy numbers.
    :return: A dictionary mapping each observed copy number to its multiplicity.
    """
    multiplicities = Counter()

    for copy_numbers in observed_copy_numbers.values():
        multiplicities.update(copy_numbers)

    multiplicities2 = count_CN_multiplicities(observed_copy_numbers)
    assert(multiplicities==multiplicities2)

    return multiplicities

def test_count_CN_multiplicities():
    observed_CNs = {
        "type1": [1, 2, 2, 3, 3, 3],
        "type2": [2, 2, 3, 3, 4],
    }

    expected_multiplicities = {
        1: 1,
        2: 4,
        3: 5,
        4: 1,
    }

    observed_multiplicities = count_CN_multiplicities(observed_CNs)

    assert observed_multiplicities == expected_multiplicities, f"Expected {expected_multiplicities}, but got {observed_multiplicities}"


def test_count_copy_number_multiplicities():
    observed_copy_numbers = {
        "type1": [1, 2, 2, 3, 3, 3],
        "type2": [2, 2, 3, 3, 4],
    }

    expected_multiplicities = {
        1: 1,
        2: 4,
        3: 5,
        4: 1,
    }

    observed_multiplicities = count_copy_number_multiplicities(observed_copy_numbers)

    assert observed_multiplicities == expected_multiplicities, f"Expected {expected_multiplicities}, but got {observed_multiplicities}"


test_count_CN_multiplicities()
test_count_copy_number_multiplicities()




def random_decimal(min_value, max_value, decimal_places):
    """
    Generate a random decimal number between min_value and max_value with a specified number of decimal places.
    
    Args:
    min_value (float): The minimum value for the generated number.
    max_value (float): The maximum value for the generated number.
    decimal_places (int): The number of decimal places for the generated number.
    
    Returns:
    float: A random decimal number between min_value and max_value with the specified number of decimal places.
    """
    random_number = random.uniform(min_value, max_value)
    return round(random_number, decimal_places)


def random_integer_log_scale(min_value, max_value):
    """
    Generate a random integer between min_value and max_value on a logarithmic scale.
    
    Args:
    min_value (int): The minimum value for the generated number.
    max_value (int): The maximum value for the generated number.
    
    Returns:
    int: A random integer between min_value and max_value on a logarithmic scale.
    """
    log_min = math.log(min_value)
    log_max = math.log(max_value)
    random_log = random.uniform(log_min, log_max)
    return int(round(math.exp(random_log)))


def biased_sample(p, min_value, max_value):
    """
    Generate a biased random integer between min_value and max_value, with a specified probability of returning the min_value.
    
    Args:
    p (float): The probability of returning the min_value.
    min_value (int): The minimum value for the generated number.
    max_value (int): The maximum value for the generated number.
    
    Returns:
    int: A biased random integer between min_value and max_value.
    """
    if random.random() < p:
        return min_value
    else:
        return random.randint(min_value + 1, max_value)


def print_simulation_parameters(pre, mid, post, p_up, p_down, rate):
    logging.info("SIMULATION PARAMETERS ARE:")
    logging.info("pre: %s", pre)
    logging.info("mid: %s", mid)
    logging.info("post: %s", post)
    logging.info("p_up: %s", p_up)
    logging.info("p_down: %s", p_down)
    logging.info("rate: %s", rate)
    logging.info("")


def print_simulated_genome_data(simulated_chromosomes, truth_trees, CN_trees, observed_CNs, observed_CN_multiplicities):
    logging.info("Simulated genome was:")
    for chrom in simulated_chromosomes:
        logging.debug("chrom: %s", str(chrom))
        logging.debug("%s", simulated_chromosomes[chrom])

    logging.debug("the truth trees are:")
    for chrom_type in truth_trees:
        logging.debug("%s", truth_trees[chrom_type])

    
    for chrom_type in CN_trees:
        logging.info("the CN simplified trees are:")
        logging.info("%s", CN_trees[chrom_type])

        logging.info("observed chromosomal copynumbers")
        logging.info("%s", observed_CNs[chrom_type])


def generate_dirichlet_probability(alpha):
    """
    Generate probabilities based on the Dirichlet distribution.

    Args:
        alpha: list or np.array
            Parameters of the Dirichlet distribution.

    Returns:
        np.array: Probabilities generated from the Dirichlet distribution, 
                  rounded to 2 decimal places and adjusted to sum to 1.
    """
    probabilities = np.random.dirichlet(alpha)

    # Round probabilities to 2 decimal places
    probabilities = np.round(probabilities, 2)

    # Adjust the last probability so the sum is 1
    probabilities[-1] = 1 - np.sum(probabilities[:-1])

    return probabilities


def test_generate_dirichlet_probability():
    alpha = [20.0, 20.0, 100.0] # [0.5, 0.5, 0.5]
    probabilities = generate_dirichlet_probability(alpha)

    assert isinstance(probabilities, np.ndarray), f"Output is not a numpy array: {type(probabilities)}"
    assert np.isclose(probabilities.sum(), 1, atol=1e-2), f"Probabilities do not sum to 1 within tolerance: {probabilities.sum()}"
    assert (probabilities >= 0).all() and (probabilities <= 1).all(), f"Probabilities contain values outside [0, 1]: {probabilities}"
    assert np.isclose(np.round(probabilities, 2).sum(), 1, atol=1e-8), "Rounded probabilities do not sum to 1"
    #assert np.round(probabilities, 2).sum() == 1, "Rounded probabilities do not sum to 1"


test_generate_dirichlet_probability()

def generate_poisson(min_value, max_value, lam):
    while True:
        value = np.random.poisson(lam)
        if min_value <= value <= max_value:
            return value

def assign_values(p1, p2, p3, min_value, max_value, lam):
    """
    Assign -1 to mid and post, post only, or neither, based on the provided probabilities.

    Args:
    p1 (float): The probability of mid and post both being -1.
    p2 (float): The probability of post being -1 and mid being not -1.
    p3 (float): The probability of mid and post both being not -1.

    min_value (int): The minimum value for the generated number.
    max_value (int): The maximum value for the generated number.
    lam (float): The rate parameter for the Poisson distribution.

    Returns:
    tuple of int: A tuple of two integers.
    """
    choices = [(generate_poisson(min_value, max_value, lam), -1),
               (-1, -1),
               (generate_poisson(min_value, max_value, lam), generate_poisson(min_value, max_value, lam))]

    return random.choices(choices, weights=[p2, p1, p3], k=1)[0]


def simulate_parameters_not_given_as_arguments(args):
    if args.max_epochs is None:
        args.max_epochs = 8

    max_anue = (args.max_epochs-2)//3

    if args.lam is None:
        args.lam = 1  # Default value for the rate parameter of the Poisson distribution

    if args.pre is None:
        args.pre = generate_poisson(0, max_anue, args.lam)

    if args.mid is None or args.post is None:
        args.mid, args.post = assign_values(0.4, 0.4, 0.2, -1, max_anue, args.lam)

    args.total_epochs = args.pre + args.mid + args.post + 2

    if args.p_up is None or args.p_down is None:
        if args.alpha is None:
            args.alpha = [2, 2, 10]  # Default values
        probabilities = generate_dirichlet_probability(args.alpha)
        args.p_up, args.p_down, _ = probabilities

    if args.rate is None:
        args.rate = random_integer_log_scale(10000,1000000)

    print_simulation_parameters(
        pre=args.pre,
        mid=args.mid,
        post=args.post,
        p_up=args.p_up,
        p_down=args.p_down,
        rate=args.rate
    )   

    # Include the simulation_filename in the return statement
    return args

import numpy as np
from scipy.stats import poisson
from typing import Tuple, Optional, List
from dataclasses import dataclass


@dataclass
class SimulationArgs:
    """Class to encapsulate simulation arguments."""
    max_epochs: Optional[int] = None
    lam: Optional[float] = None
    pre: Optional[int] = None
    mid: Optional[int] = None
    post: Optional[int] = None
    p_up: Optional[float] = None
    p_down: Optional[float] = None
    rate: Optional[int] = None
    alpha: Optional[List[int]] = None
    total_epochs: Optional[int] = None
    prob_gd_rounds: Optional[List[float]] = None


class EpochAssigner:
    """The class responsible for assigning values for pre, mid, and post epochs."""
    def __init__(self, num_gd_rounds: int, lam: float, max_epochs: int) -> None:
        """
        Initialize the EpochAssigner with the number of genome doubling rounds,
        Poisson parameter lambda, and max epochs.

        :param num_gd_rounds: Number of genome doubling rounds.
        :param lam: Poisson distribution parameter lambda for generating the number of epochs.
        :param max_epochs: Maximum number of epochs.
        """
        self.num_gd_rounds = num_gd_rounds
        self.prob_E = self._compute_prob_E(lam, max_epochs)
        self.pre, self.mid, self.post = 0, 0, 0
        self.max_epochs = max_epochs

    @staticmethod
    def _compute_prob_E(lam: float, max_epochs: int) -> np.ndarray:
        """Compute the Poisson probabilities and normalize them.

        :param lam: Poisson distribution parameter.
        :param max_epochs: Maximum number of epochs.
        :return: Normalized Poisson probabilities.
        """
        prob_E = [poisson.pmf(i, lam) for i in range(max_epochs+1)]
        prob_E = prob_E / np.sum(prob_E)
        return prob_E

    def _get_random_choice(self) -> int:
        """Returns a random choice considering the calculated probabilities."""
        return np.random.choice(range(len(self.prob_E)), p=self.prob_E)

    def assign_values(self) -> None:
        """Assign values to the pre, mid, and post epochs based on the number of genome doubling rounds."""
        if self.num_gd_rounds == 0:
            self._assign_values_for_no_rounds()
        elif self.num_gd_rounds == 1:
            self._assign_values_for_one_round()
        elif self.num_gd_rounds == 2:
            self._assign_values_for_two_rounds()

    def _get_valid_random_choice(self, number_of_periods: int) -> Tuple[int, ...]:
        """Return valid random choices that do not exceed a given limit."""
        choices = tuple(self._get_random_choice() for _ in range(number_of_periods))
        while sum(choices) >= len(self.prob_E) - number_of_periods +1:
            choices = tuple(self._get_random_choice() for _ in range(number_of_periods))
        return choices

    def _assign_values_for_no_rounds(self) -> None:
        """Assign values for the case when there are no genome doubling rounds."""
        self.pre, = self._get_valid_random_choice(number_of_periods=1)
        self.mid = self.post = -1

    def _assign_values_for_one_round(self) -> None:
        """Assign values for the case when there is one genome doubling round."""
        self.pre, self.mid = self._get_valid_random_choice(number_of_periods=2)
        self.post = -1

    def _assign_values_for_two_rounds(self) -> None:
        """Assign values for the case when there are two genome doubling rounds."""
        self.pre, self.mid, self.post = self._get_valid_random_choice(number_of_periods=3)


def simulate_parameters_not_given_as_arguments(args: SimulationArgs) -> SimulationArgs:
    """Simulates default parameters if not provided in the arguments.

    :param args: Simulation arguments object.
    :return: Updated simulation arguments with default values filled in if not provided.
    """
    args.max_epochs = args.max_epochs if args.max_epochs is not None else 8

    args.lam = args.lam if args.lam is not None else 1  # Default value for the rate parameter of the Poisson distribution

    num_gd_rounds = np.random.choice([0, 1, 2], p=args.prob_gd_rounds) if args.prob_gd_rounds is not None else np.random.choice([0, 1, 2], p=[0.4, 0.4, 0.2])

    if args.mid is None or args.post is None:
        epoch_assigner = EpochAssigner(num_gd_rounds, args.lam, args.max_epochs)
        epoch_assigner.assign_values()
        args.pre, args.mid, args.post = epoch_assigner.pre, epoch_assigner.mid, epoch_assigner.post

    args.total_epochs = args.pre + args.mid + args.post + 2

    if args.p_up is None or args.p_down is None:
        if args.alpha is None:
            args.alpha = [2, 2, 10]  # Default values
        probabilities = generate_dirichlet_probability(args.alpha)
        args.p_up, args.p_down, _ = probabilities

    if args.rate is None:
        args.rate = random_integer_log_scale(10000, 1000000)

    print_simulation_parameters(
        pre=args.pre,
        mid=args.mid,
        post=args.post,
        p_up=args.p_up,
        p_down=args.p_down,
        rate=args.rate
    )

    return args


def save_results_to_file(test_case, simulation_filename, simulated_chromosomes, truth_trees, CN_trees, observed_CNs, observed_CN_multiplicities, pre, mid, post, p_up, p_down, rate):
    file_name = f'{SIMULATIONS_FILE_FOLDER}/{simulation_filename}_{test_case}.pickle'
    with open(file_name, 'wb') as f:
        pickle.dump({
            'simulated_chromosomes': simulated_chromosomes,
            'truth_trees': truth_trees,
            'CN_trees': CN_trees,
            'observed_CNs': observed_CNs,
            'observed_CN_multiplicities': observed_CN_multiplicities,
            'pre': pre,
            'mid': mid,
            'post': post,
            'p_up': p_up,
            'p_down': p_down,
            'rate': rate
        }, f)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Simulate single chromosome data with Poisson timestamps."
    )

    parser.add_argument(
        "-t", 
        "--test_case", 
        type=int, 
        default=0, 
        help="Name of the test case you want to save the results to."
    )

    parser.add_argument(
        "-f", 
        "--simulation_filename", 
        type=str, 
        default="simulation", 
        help="Name of the simulation (optional)."
    )

    parser.add_argument(
        "-b", 
        "--debug", 
        action='store_true', 
        help="Enable debug logging."
    )

    parser.add_argument(
        "-p", 
        "--pre", 
        type=int, 
        help="Pre parameter value (optional, default: random.randint(0, (max_epochs-2)//3))."
    )

    parser.add_argument(
        "-m", 
        "--mid", 
        type=int, 
        help="Mid parameter value (optional, default: biased_sample(0.5, -1, (max_epochs-2)//3))."
    )

    parser.add_argument(
        "-o", 
        "--post", 
        type=int, 
        help="Post parameter value (optional, default: biased_sample(0.8, -1, (max_epochs-2)//3))."
    )

    parser.add_argument(
        "-u", 
        "--p_up", 
        type=float, 
        help="P_up parameter value (optional, default: random_decimal(0.10, 0.30, 2))."
    )

    parser.add_argument(
        "-d", 
        "--p_down", 
        type=float, 
        help="P_down parameter value (optional, default: random_decimal(p_up * 0.5, min(0.3, p_up * 2), 2))."
    )

    parser.add_argument(
        "-r", 
        "--rate", 
        type=int, 
        help="Rate parameter value (optional, default: random integer between 100 and 100000)."
    )

    parser.add_argument(
        "-e", 
        "--max_epochs", 
        type=int, 
        default=8, 
        help="Max epochs parameter value (optional, default is 8)."
    )

    parser.add_argument(
        "-l", 
        "--lam", 
        type=float, 
        default=1.0, 
        help="Lam parameter for Dirichlet distribution (optional, default is 1.0)."
    )

    parser.add_argument(
        "-a", 
        "--alpha", 
        nargs='*', 
        type=float, 
        default=[20.0, 20.0, 100.0], 
        help="Alpha parameter for Dirichlet distribution (optional, default is [20.0, 20.0, 100.0])."
    )
    
    parser.add_argument(
        "-g",
        "--prob_gd_rounds",
        type=float,
        nargs='+',
        help="Probability of gd rounds (optional, default: None)."
    )


    return parser.parse_args()


def main():
    # get simulation parameters:
    args = parse_arguments()
    
    if args.debug:
        logging.basicConfig(level = logging.DEBUG)
    else:
        logging.basicConfig(level = logging.INFO)
        
    args = simulate_parameters_not_given_as_arguments(args)    
    
    # simulate the genome
    simulated_chromosomes = simulate_cancer_genome(
        p_up = args.p_up,
        p_down = args.p_down,
        pre = args.pre,
        mid = args.mid,
        post = args.post,
        rate = args.rate
    )   
    
    # do some basic analysis of the simulated genome to describe it
    truth_trees = create_truth_trees_from_simulation(
        simulated_chromosomes = simulated_chromosomes,
        max_epochs = args.pre + args.mid + args.post + 2
    )   
    CN_trees = convert_truth_trees_to_CN_trees(truth_trees = copy.deepcopy(truth_trees))
    observed_CNs = count_copy_numbers(simulated_chromosomes = simulated_chromosomes)
    observed_CN_multiplicities = count_copy_number_multiplicities(observed_copy_numbers = observed_CNs)

    print_simulated_genome_data(
        simulated_chromosomes = simulated_chromosomes,
        truth_trees = truth_trees,
        CN_trees = CN_trees,
        observed_CNs = observed_CNs,
        observed_CN_multiplicities = observed_CN_multiplicities
    )

    # save the results
    save_results_to_file(
        test_case = args.test_case,
        simulation_filename = args.simulation_filename,
        simulated_chromosomes = simulated_chromosomes,
        truth_trees = truth_trees,
        CN_trees = CN_trees,
        observed_CNs = observed_CNs,
        observed_CN_multiplicities = observed_CN_multiplicities,
        pre = args.pre,
        mid = args.mid,
        post = args.post,
        p_up = args.p_up,
        p_down = args.p_down,
        rate = args.rate
    )

if __name__ == "__main__":
    main()

