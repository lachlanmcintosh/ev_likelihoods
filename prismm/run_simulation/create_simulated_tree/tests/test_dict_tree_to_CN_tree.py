import pytest
from dict_tree_to_CN_tree import get_child_and_complement_trees_from_truth_tree

@pytest.fixture
def truth_tree_none():
    return {'unique_identifier': 0, 'parent': -1, 'epoch_created': 0, 'paternal': True, 'dead': False, 'child': None, 'complement': None, 'copy_number': 1, 'SNV_multiplicity': 0, 'epoch_killed': 0}

@pytest.fixture
def truth_tree_both():
    return {'unique_identifier': -1, 'parent': None, 'epoch_created': 0, 'paternal': None, 'child': {'unique_identifier': 0, 'parent': -1, 'epoch_created': 0, 'paternal': True, 'dead': False, 'child': None, 'complement': None, 'copy_number': 1, 'SNV_multiplicity': 0, 'epoch_killed': 0}, 'complement': {'unique_identifier': 23, 'parent': -1, 'epoch_created': 0, 'paternal': False, 'dead': False, 'child': None, 'complement': None, 'copy_number': 1, 'SNV_multiplicity': 0, 'epoch_killed': 0}, 'copy_number': 2, 'SNV_multiplicity': 0, 'epoch_killed': 0}

def test_get_child_and_complement_trees_from_truth_tree_none(truth_tree_none):
    assert get_child_and_complement_trees_from_truth_tree(truth_tree_none) == (None, None)

def test_get_child_and_complement_trees_from_truth_tree_both(truth_tree_both):
    assert get_child_and_complement_trees_from_truth_tree(truth_tree_both) == ([1], [1])

#################
from dict_tree_to_CN_tree import convert_truth_tree_to_CN_tree

def test_convert_truth_tree_to_CN_tree_none(truth_tree_none):
    assert convert_truth_tree_to_CN_tree(truth_tree_none) == [1]

def test_convert_truth_tree_to_CN_tree_both(truth_tree_both):
    assert convert_truth_tree_to_CN_tree(truth_tree_both) == [2, [1], [1]]

# Additional test cases for failure or edge cases can be added similarly.

import pytest
from module import convert_truth_trees_to_CN_trees

def test_convert_truth_trees_to_CN_trees_typical_case():
    """
    Test that the function works for a typical input.
    """
    truth_trees = {'type1': {'copy_number': 1, 'child': {'copy_number': 2,
        'child': None, 'complement': None}, 'complement': {'copy_number': 3,
        'child': None, 'complement': None}}}
    expected_CN_trees = {'type1': [1, [3], [2]]}
    result_CN_trees = convert_truth_trees_to_CN_trees(truth_trees)
    assert result_CN_trees == expected_CN_trees, f'Expected {expected_CN_trees}, but got {result_CN_trees}'

def test_convert_truth_trees_to_CN_trees_large_input():
    """
    Test that the function can handle a large number of trees.
    """
    truth_trees = {i: [2, [1], [1]] for i in range(23)}
    expected_CN_trees = {i: [2, [1], [1]] for i in range(23)}
    result_CN_trees = convert_truth_trees_to_CN_trees(truth_trees)
    assert result_CN_trees == expected_CN_trees, f'Expected {expected_CN_trees}, but got {result_CN_trees}'

def test_convert_truth_trees_to_CN_trees_nested_trees():
    """
    Test that the function can handle nested trees.
    """
    truth_trees = {0: [2, [2, [1], [1]], [0]], 11: [3, [2, [1], [1]], [1]], 20: [3, [2, [1], [1]], [1]]}
    expected_CN_trees = {0: [2, [2, [1], [1]], [0]], 11: [3, [2, [1], [1]], [1]], 20: [3, [2, [1], [1]], [1]]}
    result_CN_trees = convert_truth_trees_to_CN_trees(truth_trees)
    assert result_CN_trees == expected_CN_trees, f'Expected {expected_CN_trees}, but got {result_CN_trees}'

@pytest.mark.parametrize("truth_trees, expected_CN_trees", [
    ({}, {}),
    ({'type1': None}, {'type1': None}),
])
def test_convert_truth_trees_to_CN_trees_edge_cases(truth_trees, expected_CN_trees):
    """
    Test that the function can handle edge cases like an empty dictionary or None values.
    """
    result_CN_trees = convert_truth_trees_to_CN_trees(truth_trees)
    assert result_CN_trees == expected_CN_trees, f'Expected {expected_CN_trees}, but got {result_CN_trees}'

@pytest.mark.parametrize("truth_trees", [
    "invalid",
    123,
])
def test_convert_truth_trees_to_CN_trees_failures(truth_trees):
    """
    Test that the function raises an exception for invalid inputs.
    """
    with pytest.raises(TypeError):
        convert_truth_trees_to_CN_trees(truth_trees)
