 import pytest
from check_trees import are_all_descendants_zero_in_tree_structure

def test_are_all_descendants_zero_in_tree_structure_none_node():
    """
    Test when node is None.
    """
    result = are_all_descendants_zero_in_tree_structure(None)
    assert result is True, 'Expected True, but got False.'


def test_are_all_descendants_zero_in_tree_structure_zero_copy_number_no_descendants():
    """
    Test when copy_number is zero and there are no descendants.
    """
    node = {'copy_number': 0, 'child': None, 'complement': None}
    result = are_all_descendants_zero_in_tree_structure(node)
    assert result is True, 'Expected True, but got False.'


def test_are_all_descendants_zero_in_tree_structure_non_zero_copy_number_no_descendants():
    """
    Test when copy_number is non-zero and there are no descendants.
    """
    node = {'copy_number': 1, 'child': None, 'complement': None}
    result = are_all_descendants_zero_in_tree_structure(node)
    assert result is False, 'Expected False, but got True.'


def test_are_all_descendants_zero_in_tree_structure_zero_copy_number_with_descendants():
    """
    Test when copy_number is zero and there are descendants.
    """
    node = {
        'copy_number': 0, 
        'child': {'copy_number': 0, 'child': None, 'complement': None}, 
        'complement': {'copy_number': 0, 'child': None, 'complement': None}
    }
    result = are_all_descendants_zero_in_tree_structure(node)
    assert result is True, 'Expected True, but got False.'


def test_are_all_descendants_zero_in_tree_structure_non_zero_copy_number_with_descendants():
    """
    Test when copy_number is non-zero and there are descendants.
    """
    node = {
        'copy_number': 1, 
        'child': {'copy_number': 0, 'child': None, 'complement': None}, 
        'complement': {'copy_number': 0, 'child': None, 'complement': None}
    }
    result = are_all_descendants_zero_in_tree_structure(node)
    assert result is False, 'Expected False, but got True.'

#####################

from check_trees import check_and_remove_redundant_node_in_tree_structure

def test_check_and_remove_redundant_node_no_main_child():
    """
    Test when there is no main child.
    """
    node = {'unique_identifier': 0, 'parent': -1, 'epoch_created': 0, 'paternal': True, 'dead': False, 'child': None, 'complement': None, 'copy_number': 1, 'SNV_multiplicity': 0}
    result = check_and_remove_redundant_node_in_tree_structure(node, None, None)
    assert result == node, 'Expected node to remain unchanged when there is no main child.'

def test_check_and_remove_redundant_node_main_child_different_copy_number():
    """
    Test when the main child has a different copy number than the node.
    """
    node = {'unique_identifier': 0, 'parent': -1, 'epoch_created': 0, 'paternal': True, 'dead': False, 'child': None, 'complement': None, 'copy_number': 1, 'SNV_multiplicity': 0}
    main_child = {'unique_identifier': 1, 'parent': 0, 'epoch_created': 1, 'paternal': True, 'dead': False, 'child': None, 'complement': None, 'copy_number': 2, 'SNV_multiplicity': 0}
    result = check_and_remove_redundant_node_in_tree_structure(node, main_child, None)
    assert result == node, 'Expected node to remain unchanged when main child has different copy number.'

def test_check_and_remove_redundant_node_main_child_same_copy_number_no_other_child():
    """
    Test when the main child has the same copy number as the node, and there is no other child.
    """
    node = {'unique_identifier': 0, 'parent': -1, 'epoch_created': 0, 'paternal': True, 'dead': False, 'child': None, 'complement': None, 'copy_number': 1, 'SNV_multiplicity': 0}
    main_child = {'unique_identifier': 1, 'parent': 0, 'epoch_created': 1, 'paternal': True, 'dead': False, 'child': None, 'complement': None, 'copy_number': 1, 'SNV_multiplicity': 0}
    result = check_and_remove_redundant_node_in_tree_structure(node, main_child, None)
    assert result == main_child and result['epoch_created'] == node['epoch_created'], 'Expected node to be replaced with main child and have the same epoch_created when main child has the same copy number and there is no other child.'

def test_check_and_remove_redundant_node_main_child_same_copy_number_other_child_non_zero():
    """
    Test when the main child has the same copy number as the node, and the other child has a non-zero copy number.
    """
    node = {'unique_identifier': 0, 'parent': -1, 'epoch_created': 0, 'paternal': True, 'dead': False, 'child': None, 'complement': None, 'copy_number': 1, 'SNV_multiplicity': 0}
    main_child = {'unique_identifier': 1, 'parent': 0, 'epoch_created': 1, 'paternal': True, 'dead': False, 'child': None, 'complement': None, 'copy_number': 1, 'SNV_multiplicity': 0}
    other_child = {'unique_identifier': 2, 'parent': 0, 'epoch_created': 2, 'paternal': False, 'dead': False, 'child': None, 'complement': None, 'copy_number': 1, 'SNV_multiplicity': 0}
    with pytest.raises(ValueError, match="Invalid tree: child and complement copy_number do not add up to parent's copy_number"):
        check_and_remove_redundant_node_in_tree_structure(node, main_child, other_child)
