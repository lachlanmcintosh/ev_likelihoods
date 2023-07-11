import pytest
from simulate_anueploidy_agnostic import (
    simulate_anueploidy_agnostic,
    simulate_agnostic_chromosome_loss,
    simulate_agnostic_chromosome_gain
)

from initialize_simulated_chromosomes import initialize_simulated_chromosomes
import random

# Fixing the seed for repeatability
random.seed(1)

def test_input_type_failure(method, *args):
    """
    Abstracts the common logic for testing failure cases due to type errors
    """
    for i, arg in enumerate(args):
        args_list = list(args)
        args_list[i] = "incorrect_type"
        with pytest.raises(TypeError):
            method(*args_list)

def test_simulate_anueploidy_agnostic():
    chromosomes = initialize_simulated_chromosomes()
    epoch = 1
    chrom_counts = [46, 0]

    for chrom_count in chrom_counts:
        updated_chrom_count, updated_chromosomes = simulate_anueploidy_agnostic(chromosomes, epoch, chrom_count)

        if chrom_count != 0:
            assert 0 <= updated_chrom_count <= 2 * chrom_count, 'Chromosome count should be within the range [0, 2 * chrom_count] after 1 round of simulation'
            for chrom_list in updated_chromosomes.values():
                for chrom in chrom_list:
                    if chrom['epoch_created'] > epoch - 1:
                        assert chrom['epoch_created'] == epoch, 'Newly created chromosomes should have epoch_created = 1'
        else:
            assert updated_chrom_count == 0, 'With 0 initial chromosomes, updated chromosome count should be 0 as well.'

    test_input_type_failure(simulate_anueploidy_agnostic, chromosomes, epoch, 46)

def test_simulate_agnostic_chromosome_loss():
    chromosomes = [{"dead": False} for _ in range(10)]
    chrom_counts = [10, 0]

    for chrom_count in chrom_counts:
        updated_chrom_count = simulate_agnostic_chromosome_loss(chromosomes, chrom_count)

        if chrom_count != 0:
            assert any(chrom['dead'] for chrom in chromosomes), "No chromosome has been marked as dead"
            assert updated_chrom_count == chrom_count - sum(chrom['dead'] for chrom in chromosomes), "The updated chromosome count is not correctly updated"
        else:
            assert updated_chrom_count == 0, 'With 0 initial chromosomes, updated chromosome count should be 0 as well.'

    test_input_type_failure(simulate_agnostic_chromosome_loss, chromosomes, 10)

def test_simulate_agnostic_chromosome_gain():
    epoch = 1
    chrom_count = 10
    chromosomes = [{"dead": dead, "id": i, "epoch_created": 0} for i, dead in enumerate([False]*chrom_count + [True]*chrom_count)]

    for chroms in [chromosomes[:chrom_count], chromosomes[chrom_count:]]:
        updated_chrom_count = simulate_agnostic_chromosome_gain(chroms, epoch, chrom_count)

        if not all(chrom['dead'] for chrom in chroms):
            assert len(chroms) > chrom_count, "No new chromosome has been added"
            assert updated_chrom_count == len(chroms), "The updated chromosome count is not correctly updated"
        else:
            assert updated_chrom_count == chrom_count, 'With all initial chromosomes dead, updated chromosome count should not increase.'

    test_input_type_failure(simulate_agnostic_chromosome_gain, chromosomes, epoch, chrom_count)
