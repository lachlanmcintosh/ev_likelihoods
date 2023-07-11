import pytest
import numpy as np
from typing import Any, Dict, List, Tuple

# Test that the function handles empty inputs correctly
def test_empty_simulate_snvs():
    simulated_chromosomes = {}
    lengths = {}
    rate = 0.05
    epoch = 1
    snv_count = 0
    final_snv_count, updated_chromosomes = simulate_snvs(simulated_chromosomes, lengths, rate, epoch, snv_count)
    assert final_snv_count == 0
    assert updated_chromosomes == {}

# Test that the function correctly handles dead chromosomes
def test_dead_chromosomes():
    simulated_chromosomes = {1: [{"dead": True, "SNVs": []}]}
    lengths = {1: 1000}
    rate = 0.05
    epoch = 1
    snv_count = 0
    final_snv_count, updated_chromosomes = simulate_snvs(simulated_chromosomes, lengths, rate, epoch, snv_count)
    assert final_snv_count == 0
    assert updated_chromosomes == {1: [{"dead": True, "SNVs": []}]}

# Test the rate parameter edge case with rate = 0
def test_rate_zero():
    simulated_chromosomes = {1: [{"dead": False, "SNVs": []}]}
    lengths = {1: 1000}
    rate = 0
    epoch = 1
    snv_count = 0
    final_snv_count, updated_chromosomes = simulate_snvs(simulated_chromosomes, lengths, rate, epoch, snv_count)
    assert final_snv_count == 0
    assert len(updated_chromosomes[1][0]['SNVs']) == 0

# Test that the Poisson distribution behaves as expected for very large rates
def test_large_rate():
    simulated_chromosomes = {1: [{"dead": False, "SNVs": []}]}
    lengths = {1: 1}
    rate = 1e6  # Expected to result in a very large number of SNVs
    epoch = 1
    snv_count = 0
    np.random.seed(0)  # Set a random seed for reproducibility
    final_snv_count, updated_chromosomes = simulate_snvs(simulated_chromosomes, lengths, rate, epoch, snv_count)
    assert final_snv_count > 1e5  # Expect at least 100,000 SNVs for a rate of 1e6
    assert len(updated_chromosomes[1][0]['SNVs']) == final_snv_count
