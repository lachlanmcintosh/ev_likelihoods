import pytest
import numpy as np
from random_number_generator import random_decimal, random_integer_log_scale, generate_dirichlet_probability, generate_poisson

@pytest.mark.parametrize("min_value, max_value, decimal_places", [(1.0, 2.0, 2), (-2.0, -1.0, 2)])
def test_random_decimal(min_value, max_value, decimal_places):
    """
    Test that random_decimal generates a number within the given range
    and with the correct number of decimal places.
    """
    result = random_decimal(min_value, max_value, decimal_places)
    assert min_value <= result <= max_value
    assert len(str(result).split('.')[1]) == decimal_places

@pytest.mark.parametrize("min_value, max_value, decimal_places", [(2.0, 1.0, 2), (1.0, 2.0, -1), (0, 0, 2)])
def test_random_decimal_errors(min_value, max_value, decimal_places):
    """
    Test that random_decimal throws an exception for invalid ranges or decimal places.
    """
    with pytest.raises(ValueError):
        random_decimal(min_value, max_value, decimal_places)

@pytest.mark.parametrize("min_value, max_value", [(10000, 1000000), (-10000, -1000), (0, 0), (1000000, 10000)])
def test_random_integer_log_scale(min_value, max_value):
    """
    Test that random_integer_log_scale generates a number within the given range.
    """
    if min_value <= 0 or min_value > max_value:
        with pytest.raises(ValueError):
            random_integer_log_scale(min_value, max_value)
    else:
        result = random_integer_log_scale(min_value, max_value)
        assert min_value <= result <= max_value

def test_random_integer_log_scale_log_distribution():
    """
    Test that random_integer_log_scale generates numbers that follow a logarithmic distribution.
    """
    results = [random_integer_log_scale(10000, 1000000) for _ in range(10000)]
    assert sum(results) / len(results) < ((10000 + 1000000) / 2)

@pytest.mark.parametrize("alpha", [[20.0, 20.0, 100.0], [-20.0, 20.0, 100.0]])
def test_generate_dirichlet_probability(alpha):
    """
    Test that generate_dirichlet_probability returns a numpy array of probabilities that sum to 1.
    """
    if any(a < 0 for a in alpha):
        with pytest.raises(ValueError):
            generate_dirichlet_probability(alpha)
    else:
        probabilities = generate_dirichlet_probability(alpha)
        assert isinstance(probabilities, np.ndarray)
        assert np.isclose(probabilities.sum(), 1, atol=0.01)
        assert (probabilities >= 0).all() and (probabilities <= 1).all()

@pytest.mark.parametrize("min_value, max_value, lam", [(0, 2, 1.0), (-1, 2, 1.0), (-1, 2, -1.0), (2, 1, 1.0)])
def test_generate_poisson(min_value, max_value, lam):
    """
    Test that generate_poisson returns a value within the specified range.
    """
    if lam < 0 or min_value > max_value:
        with pytest.raises(ValueError):
            generate_poisson(min_value, max_value, lam)
    else:
        value = generate_poisson(min_value, max_value, lam)
        assert min_value <= value <= max_value
