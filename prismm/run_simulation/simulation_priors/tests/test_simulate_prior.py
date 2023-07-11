import random
import pytest
from argparse import Namespace
from simulate_prior import assign_values, generate_poisson, simulate_parameters_not_given_as_arguments

# A global constant for test reproducibility
RANDOM_SEED = 1

# Test fixture for common setup steps
@pytest.fixture(scope="module", autouse=True)
def setup_module():
    random.seed(RANDOM_SEED)

def test_assign_values_outcomes():
    """
    Test if the assign_values function produces the expected outcomes.
    """
    outcomes = [(-1, -1), (generate_poisson(-1, 2, 1.0), -1), (generate_poisson(-1, 2, 1.0), generate_poisson(-1, 2, 1.0))]
    for _ in range(100):
        result = assign_values(0.4, 0.4, 0.2, -1, 2, 1.0)
        assert result in outcomes, f"Unexpected outcome: {result}"

def test_assign_values_invalid_input_handling():
    """
    Test if the assign_values function correctly raises exceptions for invalid inputs.
    """
    with pytest.raises(TypeError):
        assign_values("0.4", 0.4, 0.2, -1, 2, 1.0)
    with pytest.raises(TypeError):
        assign_values(0.4, 0.4, 0.2, "-1", 2, 1.0)

def test_assign_values_probability_sum():
    """
    Test if the assign_values function correctly raises an exception when the sum of probabilities doesn't equal 1.
    """
    with pytest.raises(ValueError):
        assign_values(0.4, 0.4, 0.3, -1, 2, 1.0)

def test_simulate_parameters_no_args():
    """
    Test if the simulate_parameters_not_given_as_arguments function correctly assigns default values.
    """
    args = Namespace(max_epochs=None, lam=None, pre=None, mid=None, post=None, p_up=None, p_down=None, rate=None, alpha=None)
    result = simulate_parameters_not_given_as_arguments(args)
    assert result.max_epochs == 8, f"Expected 8, but got {result.max_epochs}"
    assert result.lam == 1, f"Expected 1, but got {result.lam}"
    assert isinstance(result.pre, int), f"Expected an integer, but got {type(result.pre)}"
    assert isinstance(result.mid, int), f"Expected an integer, but got {type(result.mid)}"
    assert isinstance(result.post, int), f"Expected an integer, but got {type(result.post)}"
    assert isinstance(result.rate, int), f"Expected an integer, but got {type(result.rate)}"
    assert isinstance(result.p_up, float), f"Expected a float, but got {type(result.p_up)}"
    assert isinstance(result.p_down, float), f"Expected a float, but got {type(result.p_down)}"

def test_simulate_parameters_with_args():
    """
    Test if the simulate_parameters_not_given_as_arguments function correctly keeps the assigned values.
    """
    args = Namespace(max_epochs=10, lam=2, pre=2, mid=None, post=None, p_up=None, p_down=None, rate=None, alpha=None)
    result = simulate_parameters_not_given_as_arguments(args)
    assert result.max_epochs == 10, f"Expected 10, but got {result.max_epochs}"
    assert result.lam == 2, f"Expected 2, but got {result.lam}"
    assert result.pre == 2, f"Expected 2, but got {result.pre}"
