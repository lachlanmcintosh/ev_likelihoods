import pytest
from typing import Dict, Optional
from add_metadata import add_copynumber_to_tree_structure, assign_epoch_killed_to_tree_structure, add_SNV_multiplicity_to_tree_structure


def test_add_copynumber_to_tree_structure_empty():
    """
    Test for an empty tree.
    """
    tree = {}
    with pytest.raises(KeyError):
        add_copynumber_to_tree_structure(tree)


def test_add_copynumber_to_tree_structure_no_child_no_complement():
    """
    Test when both the child and complement are None.
    """
    tree = {'child': None, 'complement': None, 'dead': False}
    expected_tree = {'child': None, 'complement': None, 'dead': False, 'copy_number': 1}
    result_tree = add_copynumber_to_tree_structure(tree)
    assert result_tree == expected_tree, f'Expected {expected_tree}, but got {result_tree}'


def test_add_copynumber_to_tree_structure_one_level():
    """
    Test when tree has one level of child and complement.
    """
    tree = {
        'child': {'child': None, 'complement': None, 'dead': False},
        'complement': {'child': None, 'complement': None, 'dead': True}
    }
    expected_tree = {
        'child': {'child': None, 'complement': None, 'dead': False, 'copy_number': 1},
        'complement': {'child': None, 'complement': None, 'dead': True, 'copy_number': 0},
        'copy_number': 1
    }
    result_tree = add_copynumber_to_tree_structure(tree)
    assert result_tree == expected_tree, f'Expected {expected_tree}, but got {result_tree}'


def test_add_copynumber_to_tree_structure_multi_level():
    """
    Test when tree has multiple levels of child and complement.
    """
    tree = {
        'child': {
            'child': {'child': None, 'complement': None, 'dead': False},
            'complement': {'child': None, 'complement': None, 'dead': True}
        },
        'complement': {'child': None, 'complement': None, 'dead': False}
    }
    expected_tree = {
        'child': {
            'child': {'child': None, 'complement': None, 'dead': False, 'copy_number': 1},
            'complement': {'child': None, 'complement': None, 'dead': True, 'copy_number': 0},
            'copy_number': 1
        },
        'complement': {'child': None, 'complement': None, 'dead': False, 'copy_number': 1},
        'copy_number': 2
    }
    result_tree = add_copynumber_to_tree_structure(tree)
    assert result_tree == expected_tree, f'Expected {expected_tree}, but got {result_tree}'


def test_add_copynumber_to_tree_structure_existing_copy_number():
    """
    Test when copy_number already exists in the input tree.
    """
    tree = {
        'copy_number': 1,
        'child': None,
        'complement': None,
        'dead': False
    }
    expected_tree = {
        'copy_number': 1,
        'child': None,
        'complement': None,
        'dead': False
    }
    result_tree = add_copynumber_to_tree_structure(tree)
    assert result_tree == expected_tree, f'Expected {expected_tree}, but got {result_tree}'


##############

def test_add_SNV_multiplicity_to_tree_structure_empty():
    """
    Test when the tree is empty.
    """
    with pytest.raises(KeyError):
        add_SNV_multiplicity_to_tree_structure({})


def test_add_SNV_multiplicity_to_tree_structure_no_child_no_complement():
    """
    Test when both the child and complement are None.
    """
    tree = {'copy_number': 1, 'child': None, 'complement': None, 'epoch_created': 0, 'SNVs': [{'epoch_created': 0}, {'epoch_created': 1}]}
    expected_tree = tree.copy()
    expected_tree['SNV_multiplicity'] = 2
    result_tree = add_SNV_multiplicity_to_tree_structure(tree)
    assert result_tree == expected_tree, f'Expected {expected_tree}, but got {result_tree}'


def test_add_SNV_multiplicity_to_tree_structure_with_child_complement():
    """
    Test when tree has both child and complement.
    """
    tree = {
        'copy_number': 2, 'epoch_created': 0, 'SNVs': [{'epoch_created': 0}], 
        'child': {
            'copy_number': 1, 'epoch_created': 1, 'SNVs': [{'epoch_created': 1}], 
            'child': None, 'complement': None
        },
        'complement': {
            'copy_number': 1, 'epoch_created': 1, 'SNVs': [{'epoch_created': 2}], 
            'child': None, 'complement': None
        }
    }
    expected_tree = {
        'copy_number': 2, 'epoch_created': 0, 'SNV_multiplicity': 1, 
        'SNVs': [{'epoch_created': 0}], 
        'child': {
            'copy_number': 1, 'epoch_created': 1, 'SNV_multiplicity': 1, 
            'SNVs': [{'epoch_created': 1}], 'child': None, 'complement': None
        },
        'complement': {
            'copy_number': 1, 'epoch_created': 1, 'SNV_multiplicity': 1, 
            'SNVs': [{'epoch_created': 2}], 'child': None, 'complement': None
        }
    }
    result_tree = add_SNV_multiplicity_to_tree_structure(tree)
    assert result_tree == expected_tree, f'Expected {expected_tree}, but got {result_tree}'


def test_add_SNV_multiplicity_to_tree_structure_zero_copy_number():
    """
    Test when copy_number of tree is zero.
    """
    tree = {'copy_number': 0, 'child': None, 'complement': None, 'SNVs': [], 'epoch_created': 0}
    expected_tree = {'copy_number': 0, 'child': None, 'complement': None, 'SNVs': [], 'epoch_created': 0, 'SNV_multiplicity': None}
    result_tree = add_SNV_multiplicity_to_tree_structure(tree)
    assert result_tree == expected_tree, f'Expected {expected_tree}, but got {result_tree}'


def test_add_SNV_multiplicity_to_tree_structure_existing_SNV_multiplicity():
    """
    Test when SNV_multiplicity already exists in the input tree.
    """
    tree = {'copy_number': 1, 'child': None, 'complement': None, 'epoch_created': 0, 'SNVs': [{'epoch_created': 0}, {'epoch_created': 1}], 'SNV_multiplicity': 1}
    expected_tree = tree.copy()
    expected_tree['SNV_multiplicity'] = 2
    result_tree = add_SNV_multiplicity_to_tree_structure(tree)
    assert result_tree == expected_tree, f'Expected {expected_tree}, but got {result_tree}'


##############

def test_assign_epoch_killed_to_tree_structure_existing():
    """
    Test when the tree already has 'epoch_killed' values.
    """
    tree = {
        'epoch_created': 1, 
        'child': {
            'epoch_created': 2, 
            'child': None, 
            'complement': None
        }, 
        'complement': {
            'epoch_created': 2, 
            'child': None, 
            'complement': None
        }
    }
    max_epochs = 3
    expected_tree = {
        'epoch_created': 1, 
        'epoch_killed': 2, 
        'child': {
            'epoch_created': 2, 
            'epoch_killed': 3, 
            'child': None, 
            'complement': None
        }, 
        'complement': {
            'epoch_created': 2, 
            'epoch_killed': 3,
            'child': None, 
            'complement': None
        }
    }
    result_tree = assign_epoch_killed_to_tree_structure(tree, max_epochs)
    assert result_tree == expected_tree, f'Expected {expected_tree}, but got {result_tree}'


def test_assign_epoch_killed_to_tree_structure_empty():
    """
    Test when the tree is None.
    """
    result_tree = assign_epoch_killed_to_tree_structure(None, 5)
    assert result_tree is None, 'Expected None, but got a value.'


def test_assign_epoch_killed_to_tree_structure_no_child_complement():
    """
    Test when both the child and complement are None.
    """
    tree = {'epoch_created': 1, 'child': None, 'complement': None}
    expected_tree = {'epoch_created': 1, 'child': None, 'complement': None, 'epoch_killed': 3}
    result_tree = assign_epoch_killed_to_tree_structure(tree, 3)
    assert result_tree == expected_tree, f'Expected {expected_tree}, but got {result_tree}'


def test_assign_epoch_killed_to_tree_structure_child_complement_mismatch():
    """
    Test when the 'epoch_created' values of the child and complement do not match.
    """
    tree = {
        'epoch_created': 1, 
        'child': {
            'epoch_created': 2, 
            'child': None, 
            'complement': None
        }, 
        'complement': {
            'epoch_created': 3, 
            'child': None, 
            'complement': None
        }
    }
    with pytest.raises(ValueError):
        assign_epoch_killed_to_tree_structure(tree, 3)
