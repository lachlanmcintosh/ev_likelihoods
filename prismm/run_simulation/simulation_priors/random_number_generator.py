import random
import math
import numpy as np

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

def generate_poisson(max_value: int, lam: float, max_attempts: int = 1000) -> int:
    """
    Generate a random number from a Poisson distribution, 
    where the number is within the range defined by min_value and max_value.

    Args:
        min_value: The minimum value for the generated number.
        max_value: The maximum value for the generated number.
        lam: The lambda parameter for the Poisson distribution.
        max_attempts: The maximum number of attempts to generate a suitable number.

    Returns:
        A random number from a Poisson distribution within the specified range.

    Raises:
        ValueError: If a suitable number can't be generated within max_attempts.
    """
    for _ in range(max_attempts):
        value = np.random.poisson(lam)
        if 0 <= value <= max_value:
            return value
    raise ValueError(f"Could not generate a suitable random number within {max_attempts} attempts.")




