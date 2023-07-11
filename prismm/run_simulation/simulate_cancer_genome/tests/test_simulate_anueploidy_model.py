import pytest
import numpy as np
from simulate_anueploidy_model import simulate_anueploidy, simulate_chromosome_changes, initialize_simulated_chromosomes

np.random.seed(1)  # Ensure repeatability of the tests

@pytest.mark.parametrize("simulate_func", [simulate_anueploidy, simulate_chromosome_changes])
def test_simulation_func(simulate_func):
    # Initialize
    chromosomes = initialize_simulated_chromosomes()
    chrom_key = "autosome" if simulate_func == simulate_chromosome_changes else None
    chromosomes = chromosomes[chrom_key] if chrom_key else chromosomes
    epoch = 1
    chrom_count = 46

    # Typical case
    p_up = 0.1
    p_down = 0.1
    new_chromosomes, updated_chrom_count = simulate_func(chromosomes, epoch, chrom_count, p_up, p_down)
    
    assert updated_chrom_count >= chrom_count, 'Chromosome count should increase or stay the same after simulation'
    chrom_list = [new_chromosomes] if simulate_func == simulate_chromosome_changes else new_chromosomes.values()
    for chrom in chrom_list:
        if chrom['epoch_created'] > epoch - 1:
            assert chrom['epoch_created'] == 1, 'Newly created chromosomes should have epoch_created = 1'
            
    # Edge case: Zero probabilities
    p_up = 0
    p_down = 0
    new_chromosomes, updated_chrom_count = simulate_func(chromosomes, epoch, chrom_count, p_up, p_down)
    assert updated_chrom_count == chrom_count, 'Chromosome count should stay the same when probabilities are 0'

    # Edge case: All probabilities are 1
    p_up = 1
    p_down = 1
    with pytest.raises(ValueError):
        new_chromosomes, updated_chrom_count = simulate_func(chromosomes, epoch, chrom_count, p_up, p_down)

    # Failure case: Incorrect types
    invalid_cases = [("invalid", p_down), (p_up, "invalid"), ("invalid", p_up, p_down)]
    for case in invalid_cases:
        with pytest.raises(AssertionError):
            new_chromosomes, updated_chrom_count = simulate_func(chromosomes, epoch, chrom_count, *case)
    
    with pytest.raises(AssertionError):
        new_chromosomes, updated_chrom_count = simulate_func(chromosomes, "invalid", chrom_count, p_up, p_down)
    
    with pytest.raises(AssertionError):
        new_chromosomes, updated_chrom_count = simulate_func("invalid", epoch, chrom_count, p_up, p_down)
