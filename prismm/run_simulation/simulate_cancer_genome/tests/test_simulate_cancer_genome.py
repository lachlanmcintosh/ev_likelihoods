import pytest
from simulate_cancer_genome import simulate_cancer_genome

# Define a fixture to setup common parameters
@pytest.fixture
def default_params():
    return {
        'p_up': 0.1, 
        'p_down': 0.1,
        'pre': 2, 
        'mid': 1, 
        'post': 3, 
        'rate': 0.01,
        'agnostic': False
    }

# Test for default conditions
def test_simulate_cancer_genome_default(default_params):
    simulated_chromosomes = simulate_cancer_genome(**default_params)
    assert isinstance(simulated_chromosomes, list), 'Output should be a list'
    assert simulated_chromosomes, 'There should be chromosomes in the simulation'
    assert check_all_chrs_are_unique(simulated_chromosomes), 'All chromosomes in the simulation should be unique'
    assert check_expected_keys_in_simulated_chromosomes_present(simulated_chromosomes)

# Test for when mid=-1, post!=-1
def test_simulate_cancer_genome_invalid_mid_post(default_params):
    default_params.update({'mid': -1, 'post': 0})
    with pytest.raises(AssertionError):
        simulate_cancer_genome(**default_params)

# Test for when mid=0, post=-1
def test_simulate_cancer_genome_mid_zero_post_negative(default_params):
    default_params.update({'mid': 0, 'post': -1})
    simulated_chromosomes = simulate_cancer_genome(**default_params)
    assert simulated_chromosomes[-1] == 'G', "If a genome doubling round just occurred, the last element of the simulated_chromosomes should be 'G'"

# Test for when post=0
def test_simulate_cancer_genome_post_zero(default_params):
    default_params.update({'post': 0})
    simulated_chromosomes = simulate_cancer_genome(**default_params)
    assert simulated_chromosomes[-1] == 'G', "If a genome doubling round just occurred, the last element of the simulated_chromosomes should be 'G'"
