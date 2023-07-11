import pytest
from count_multiplicities import count_copy_numbers
from mymodule import count_copy_number_multiplicities
from collections import Counter

# Test data for count_copy_numbers
test_data_count_copy_numbers = [
    # Typical usage
    ({
        'type1': [{'paternal': True, 'dead': False}, {'paternal': True, 'dead': True}, {'paternal': False, 'dead': False}, {'paternal': False, 'dead': False}],
        'type2': [{'paternal': True, 'dead': False}, {'paternal': True, 'dead': False}, {'paternal': False, 'dead': True}, {'paternal': False, 'dead': True}]
    }, {'type1': [1, 2], 'type2': [2, 0]}),
    # Empty input
    ({}, {}),
    # No living chromosomes
    ({
        'type1': [{'paternal': True, 'dead': True}, {'paternal': False, 'dead': True}],
        'type2': [{'paternal': True, 'dead': True}, {'paternal': False, 'dead': True}]
    }, {'type1': [0, 0], 'type2': [0, 0]}),
    # No dead chromosomes
    ({
        'type1': [{'paternal': True, 'dead': False}, {'paternal': False, 'dead': False}],
        'type2': [{'paternal': True, 'dead': False}, {'paternal': False, 'dead': False}]
    }, {'type1': [1, 1], 'type2': [1, 1]}),
]

@pytest.mark.parametrize("simulated_chromosomes, expected", test_data_count_copy_numbers)
def test_count_copy_numbers(simulated_chromosomes, expected):
    observed = count_copy_numbers(simulated_chromosomes)
    assert observed == expected, f'Expected {expected}, but got {observed}'

# Test data for count_copy_number_multiplicities
test_data_count_copy_number_multiplicities = [
    # Normal case
    ({'type1': [1, 2, 2, 3, 3, 3], 'type2': [2, 2, 3, 3, 4]}, Counter({1: 1, 2: 4, 3: 5, 4: 1})),
    # Empty dict
    ({}, Counter()),
    # Single type same values
    ({'type1': [2, 2, 2, 2]}, Counter({2: 4})),
    # Non-integer values
    ({'type1': [1.5, 2.3, 2.3, 3.1, 3.1, 3.1]}, Counter({1.5: 1, 2.3: 2, 3.1: 3})),
    # Type with empty list
    ({'type1': [], 'type2': [2, 2, 3, 3, 4]}, Counter({2: 2, 3: 2, 4: 1})),
]

@pytest.mark.parametrize("observed_copy_numbers, expected", test_data_count_copy_number_multiplicities)
def test_count_copy_number_multiplicities(observed_copy_numbers, expected):
    observed = count_copy_number_multiplicities(observed_copy_numbers)
    assert observed == expected, f'Expected {expected}, but got {observed}'
