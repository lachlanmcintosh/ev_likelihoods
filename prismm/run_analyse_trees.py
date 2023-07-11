import pandas as pd
import sys
from general_functions import *
import math
from constants import *
import pickle
import os
import json
import argparse
from abc import ABC, abstractmethod
from typing import Optional, Callable, Dict, Tuple, List, Type, Union
import math
import pytest



def count_nodes(tree):
    """
    Counts the number of nodes in a given tree.
    
    :param tree: A dictionary representing a tree with keys 'child' and 'complement'.
    :return: The number of nodes in the tree.
    """
    if "child" in tree and tree["child"] is not None:
        return count_nodes(tree["child"]) + count_nodes(tree["complement"]) + 1
    return 1


def test_count_nodes():
    tree1 = {"child": None, "complement": None}
    tree2 = {"child": {"child": None, "complement": None}, "complement": {"child": None, "complement": None}}
    
    assert count_nodes(tree1) == 1, "Test Case 1 Failed"
    assert count_nodes(tree2) == 3, "Test Case 2 Failed"
    


def sort_tree_by_copynumber(tree):
    """
    Sorts a given tree by 'copy_number' to ensure that 'child' has a greater 'copy_number' than 'complement'.
    In case of ties, it ensures that 'child' has a smaller 'epoch_created'.
    
    :param tree: A dictionary representing a tree with keys 'child', 'complement', 'copy_number', and 'epoch_index'.
    :return: The sorted tree.
    """

    if tree is None:
        return None

    if 'child' in tree and tree['child'] is not None:
        sort_tree_by_copynumber(tree['child'])
        
    if 'complement' in tree and tree['complement'] is not None:
        sort_tree_by_copynumber(tree['complement'])

    if 'child' in tree and tree['child'] is not None and \
       'complement' in tree and tree['complement'] is not None:
        if tree['complement'].get('copy_number') > tree['child'].get('copy_number'):
            tree['child'], tree['complement'] = tree['complement'], tree['child']
        elif tree['complement'].get('copy_number') == tree['child'].get('copy_number') and \
             tree['child'].get('epoch_index') > tree['complement'].get('epoch_index'):
            tree['child'], tree['complement'] = tree['complement'], tree['child']
    
    return tree


def test_sort_tree_by_copynumber():
    tree1 = {
        "child": {
            "child": None,
            "complement": None,
            "copy_number": 2,
            "epoch_index": 1
        },
        "complement": {
            "child": None,
            "complement": None,
            "copy_number": 1,
            "epoch_index": 2
        }
    }

    sorted_tree1 = {
        "child": {
            "child": None,
            "complement": None,
            "copy_number": 2,
            "epoch_index": 1
        },
        "complement": {
            "child": None,
            "complement": None,
            "copy_number": 1,
            "epoch_index": 2
        }
    }

    tree2 = {
        "child": {
            "child": None,
            "complement": None,
            "copy_number": 1,
            "epoch_index": 2
        },
        "complement": {
            "child": None,
            "complement": None,
            "copy_number": 2,
            "epoch_index": 1
        }
    }

    sorted_tree2 = {
        "child": {
            "child": None,
            "complement": None,
            "copy_number": 2,
            "epoch_index": 1
        },
        "complement": {
            "child": None,
            "complement": None,
            "copy_number": 1,
            "epoch_index": 2
        }
    }

    assert sort_tree_by_copynumber(tree1) == sorted_tree1, "Test Case 1 Failed"
    assert sort_tree_by_copynumber(tree2) == sorted_tree2, "Test Case 2 Failed"


# Run test functions
test_count_nodes()
test_sort_tree_by_copynumber()

def are_trees_identical(tree1, tree2):
    """
    This function takes two trees as input, each represented as a dictionary, and returns a Boolean indicating whether they are
    topologically identical and have the same epoch_created value and copy_number value at each node.
    The function first checks if the unique_identifier, epoch_created, and copy_number values are equal between the two trees.
    If both trees have both child and complement keys missing, the function returns True.
    If only one tree has a child or complement key missing, the function returns False.
    If both trees have child or complement keys, the function recursively calls itself on the child or complement of the two trees.
    If all checks return True, the function returns True, indicating that the two trees are topologically identical and have the same
    epoch_created value and copy_number value at each node.
    """
    if tree1['epoch_index'] != tree2['epoch_index']:
        return False

    if tree1['copy_number'] != tree2['copy_number']:
        return False

    if (tree1.get('child') is None and
        tree2.get('child') is None and
        tree1.get('complement') is None and
        tree2.get('complement') is None):
        return True

    if (tree1.get('child') is None) != (tree2.get('child') is None):
        return False

    if (tree1.get('complement') is None) != (tree2.get('complement') is None):
        return False

    if (tree1.get('child') is not None and
            not are_trees_identical(tree1['child'], tree2['child'])):
        return False

    if (tree1.get('complement') is not None and
            not are_trees_identical(tree1['complement'], tree2['complement'])):
        return False

    return True


def test_are_trees_identical():
    tree1 = {'epoch_index': 1, 'copy_number': 1}
    tree2 = {'epoch_index': 1, 'copy_number': 2}
    tree3 = {'epoch_index': 2, 'copy_number': 1}

    test_cases = [
        (tree1, tree1, True),
        (tree1, tree2, False),
        (tree1, tree3, False),
    ]

    for t1, t2, expected in test_cases:
        assert are_trees_identical(t1, t2) == expected, f"Expected {expected}, but got {are_trees_identical(t1, t2)} for {t1} and {t2}"

    tree1['child'] = {'epoch_index': 2, 'copy_number': 1}
    tree2['child'] = {'epoch_index': 2, 'copy_number': 1}
    tree3['child'] = {'epoch_index': 3, 'copy_number': 1}

    test_cases = [
        (tree1, tree1, True),
        (tree1, tree2, False),
        (tree1, tree3, False),
    ]

    for t1, t2, expected in test_cases:
        assert are_trees_identical(t1, t2) == expected, f"Expected {expected}, but got {are_trees_identical(t1, t2)} for {t1} and {t2} with child"

    tree1['complement'] = {'epoch_index': 3, 'copy_number': 1}
    tree2['complement'] = {'epoch_index': 3, 'copy_number': 1}
    tree3['complement'] = {'epoch_index': 4, 'copy_number': 1}

    test_cases = [
        (tree1, tree1, True),
        (tree1, tree2, False),
        (tree2, tree3, False),
    ]

    for t1, t2, expected in test_cases:
        assert are_trees_identical(t1, t2) == expected, f"Expected {expected}, but got {are_trees_identical(t1, t2)} for {t1} and {t2} with complement"

    tree1['child']['complement'] = {'epoch_index': 4, 'copy_number': 1}
    tree2['child']['complement'] = {'epoch_index': 4, 'copy_number': 1}
    tree3['child']['complement'] = {'epoch_index': 5, 'copy_number': 1}

    test_cases = [
        (tree1, tree1, True),
        (tree1, tree2, False),
        (tree1, tree3, False),
    ]

    for t1, t2, expected in test_cases:
        assert are_trees_identical(t1, t2) == expected, f"Expected {expected}, but got {are_trees_identical(t1, t2)} for {t1} and {t2} with child and complement"

test_are_trees_identical()









def sum_tree_distance(tree1, tree2, diff_struct_is_inf=False):
    """
    Calculate the sum of absolute differences between the epoch_created values of each node for the nodes
    that have identical copy_number values in two given trees.

    Args:
        tree1 (dict): The first tree represented as a dictionary.
        tree2 (dict): The second tree represented as a dictionary.
        diff_struct_is_inf (bool, optional): If True, returns infinity if there is a structural difference
                                             between the trees. Default is False.

    Returns:
        float: The sum of absolute differences between the epoch_created values of nodes with identical
               copy_number values.
    """
    total = 0

    if (tree1 is not None and tree2 is not None and
            tree1['copy_number'] == tree2['copy_number'] and
            tree1["epoch_index"] is not None and
            tree2["epoch_index"] is not None):
        total += abs(tree1['epoch_index'] - tree2['epoch_index'])

    if tree1.get('child') is not None and tree2.get('child') is not None:
        total += sum_tree_distance(tree1['child'], tree2['child'], diff_struct_is_inf)
    elif diff_struct_is_inf:
        return float('inf')

    if tree1.get('complement') is not None and tree2.get('complement') is not None:
        total += sum_tree_distance(tree1['complement'], tree2['complement'], diff_struct_is_inf)
    elif diff_struct_is_inf:
        return float('inf')

    return total


def test_sum_tree_distance():
    def assert_with_message(value1, value2, message):
        assert value1 == value2, f"{message}: expected {value2}, but got {value1}"

    tree1 = {
        "copy_number": 1,
        "epoch_index": 1,
        "child": {
            "copy_number": 2,
            "epoch_index": 2,
        },
        "complement": {
            "copy_number": 3,
            "epoch_index": 3,
        }
    }

    tree2 = {
        "copy_number": 1,
        "epoch_index": 4,
        "child": {
            "copy_number": 2,
            "epoch_index": 5,
        },
        "complement": {
            "copy_number": 3,
            "epoch_index": 6,
        }
    }

    assert_with_message(sum_tree_distance(tree1, tree2), 9, "Test case 1 failed")
    assert_with_message(sum_tree_distance(tree1, tree2, diff_struct_is_inf=True), float('inf'), "Test case 2 failed")
    
    tree2["complement"] = None
    assert_with_message(sum_tree_distance(tree1, tree2, diff_struct_is_inf=True), float('inf'), "Test case 3 failed")
    assert_with_message(sum_tree_distance(tree1, tree2), 6, "Test case 4 failed")
    
    tree2["child"] = None
    assert_with_message(sum_tree_distance(tree1, tree2, diff_struct_is_inf=True), float('inf'), "Test case 5 failed")
    assert_with_message(sum_tree_distance(tree1, tree2), 3, "Test case 6 failed")


test_sum_tree_distance()


def count_nodes_with_same_attributes(tree1, tree2, attributes):
    """
    This function uses recursion to iterate through both trees and count the number of nodes that have
    the same value for the specified attributes. It first checks if the current nodes in both trees have
    the attribute key, and if they do and their values are equal, it increments the count by 1.

    :param tree1: A dictionary representing the first tree
    :param tree2: A dictionary representing the second tree
    :param attributes: A list of strings representing the attributes to be compared
    :return: The count of nodes with the same attribute value in both trees
    """
    count = 0
    for attribute in attributes:
        if attribute in tree1 and attribute in tree2 and tree1[attribute] == tree2[attribute]:
            count += 1
    for child_key in ['child', 'complement']:
        if child_key in tree1 and child_key in tree2 and tree1[child_key] is not None and tree2[child_key] is not None:
            count += count_nodes_with_same_attributes(tree1[child_key], tree2[child_key], attributes)
    return count


def test_count_nodes_with_same_attributes():
    tree1 = {'copy_number': 1, 'epoch_index': 2, 'child': {'copy_number': 2, 'epoch_index': 3}}
    tree2 = {'copy_number': 1, 'epoch_index': 2, 'child': {'copy_number': 2, 'epoch_index': 4}}

    result = count_nodes_with_same_attributes(tree1, tree2, ['copy_number', 'epoch_index'])
    expected = 3
    assert result == expected, f"Expected {expected}, but got {result}"

    tree1['complement'] = {'copy_number': 3, 'epoch_index': 5}
    tree2['complement'] = {'copy_number': 3, 'epoch_index': 6}

    result = count_nodes_with_same_attributes(tree1, tree2, ['copy_number', 'epoch_index'])
    expected = 4
    assert result == expected, f"Expected {expected}, but got {result}"

test_count_nodes_with_same_attributes()

# define two special functions
count_nodes_with_same_copy_number = lambda tree1, tree2: count_nodes_with_same_attributes(tree1, tree2, ['copy_number'])
count_nodes_with_same_properties = lambda tree1, tree2: count_nodes_with_same_attributes(tree1, tree2, ['epoch_index', 'copy_number'])


def are_trees_identical_by_epoch_and_copy_number(tree1: dict, tree2: dict) -> bool:
    """
    This function takes two trees as input, each represented as a dictionary, and returns a Boolean indicating whether they
    are topologically identical and have the same epoch_created value and copy_number value at each node.

    Args:
        tree1: A dictionary representing a tree.
        tree2: A dictionary representing a tree.

    Returns:
        A boolean indicating whether the two trees are topologically identical and have the same epoch_created value and
        copy_number value at each node.
    """

    if tree1['epoch_index'] != tree2['epoch_index'] or tree1['copy_number'] != tree2['copy_number']:
        return False

    tree1_child = tree1.get('child')
    tree2_child = tree2.get('child')
    tree1_complement = tree1.get('complement')
    tree2_complement = tree2.get('complement')

    if not tree1_child and not tree2_child and not tree1_complement and not tree2_complement:
        return True

    if bool(tree1_child) != bool(tree2_child) or bool(tree1_complement) != bool(tree2_complement):
        return False

    if tree1_child and not are_trees_identical_by_epoch_and_copy_number(tree1_child, tree2_child):
        return False

    if tree1_complement and not are_trees_identical_by_epoch_and_copy_number(tree1_complement, tree2_complement):
        return False

    return True


def test_are_trees_identical_by_epoch_and_copy_number():
    # Identical trees
    tree1 = {'epoch_index': 1, 'copy_number': 2, 'child': {'epoch_index': 3, 'copy_number': 4}}
    tree2 = {'epoch_index': 1, 'copy_number': 2, 'child': {'epoch_index': 3, 'copy_number': 4}}
    assert are_trees_identical_by_epoch_and_copy_number(tree1, tree2)

    # Different epoch_index
    tree1 = {'epoch_index': 1, 'copy_number': 2, 'child': {'epoch_index': 3, 'copy_number': 4}}
    tree2 = {'epoch_index': 2, 'copy_number': 2, 'child': {'epoch_index': 3, 'copy_number': 4}}
    assert not are_trees_identical_by_epoch_and_copy_number(tree1, tree2)

    # Different copy_number
    tree1 = {'epoch_index': 1, 'copy_number': 2, 'child': {'epoch_index': 3, 'copy_number': 4}}
    tree2 = {'epoch_index': 1, 'copy_number': 3, 'child': {'epoch_index': 3, 'copy_number': 4}}
    assert not are_trees_identical_by_epoch_and_copy_number(tree1, tree2)




def convert_dict_tree_to_list(tree, total_epochs=None, is_truth=False):
    """
    Convert a tree in dictionary format to a list representation for easy visual inspection.

    :param tree: Dictionary tree representation to be converted
    :param total_epochs: Total number of epochs in the tree (default: None)
    :param is_truth: Flag to determine if the tree is a "truth" tree (default: False)
    :return: List representation of the tree
    """

    if tree is None:
        return None

    copy_number = tree.get('copy_number')
    epoch_index = tree.get('epoch_index')

    if is_truth:
        if 'child' in tree and tree["child"] is not None:
            epoch_index = tree["child"]["epoch_index"]
        else:
            epoch_index = total_epochs

    child_tree = convert_dict_tree_to_list(tree.get('child'), total_epochs, is_truth)
    complement_tree = convert_dict_tree_to_list(tree.get('complement'), total_epochs, is_truth)

    if child_tree is None:
        return [(copy_number, epoch_index)]
    return [(copy_number, epoch_index), child_tree, complement_tree]


def test_convert_dict_tree_to_list():
    tree = {
        'copy_number': 1,
        'epoch_index': 0,
        'child': {
            'copy_number': 2,
            'epoch_index': 1,
            'child': None,
            'complement': None
        },
        'complement': {
            'copy_number': 3,
            'epoch_index': 2,
            'child': None,
            'complement': None
        }
    }

    expected_result = [(1, 0), [(2, 1)], [(3, 2)]]
    result = convert_dict_tree_to_list(tree)

    assert result == expected_result, f"Assertion failed: {result} != {expected_result}"


test_convert_dict_tree_to_list()


def add_epoch_created(tree, epoch_list):
    """
    Add epoch_created key to each node in the tree.
    
    :param tree: dict, the tree structure
    :param epoch_list: list of epochs
    :return: updated tree
    """
    if not tree:
        raise ValueError("Tree must be non-empty.")

    if not epoch_list:
        raise ValueError("Epoch list must be non-empty.")
        
    epoch_iterator = iter(epoch_list)
    
    def add_epoch_to_node(node):
        try:
            node['epoch_created'] = next(epoch_iterator)
        except StopIteration:
            raise ValueError(f"Epoch list is shorter than the number of nodes in the tree. "
                             f"Tree: {tree}, Epoch List: {epoch_list}")

        if 'child' in node:
            add_epoch_to_node(node['child'])

        if 'complement' in node:
            add_epoch_to_node(node['complement'])
    
    add_epoch_to_node(tree)
    # Check if we've consumed the entire epoch_list.
    try:
        next(epoch_iterator)
        raise ValueError(f"Epoch list is longer than the number of nodes in the tree. "
                         f"Tree: {tree}, Epoch List: {epoch_list}")
    except StopIteration:
        pass

    return tree


def test_add_epoch_created():
    # Test 1: Tree with three nodes
    tree1 = {'copy_number': 2, 'child': {'copy_number': 1}, 'complement': {'copy_number': 1}}
    epoch_list1 = [1, 2, 3]
    expected_output1 = {
        'copy_number': 2,
        'epoch_created': 1,
        'child': {
            'copy_number': 1,
            'epoch_created': 2
        },
        'complement': {
            'copy_number': 1,
            'epoch_created': 3
        }
    }
    assert add_epoch_created(tree1, epoch_list1) == expected_output1, f"Error: {add_epoch_created(tree1, epoch_list1)}"

    # Test 2: More complex tree
    tree2 = {'copy_number': 6, 'child': {'copy_number': 4, 'child': {'copy_number': 3}, 'complement': {'copy_number': 1}}, 'complement': {'copy_number': 2, 'child': {'copy_number': 1}, 'complement': {'copy_number': 1}}}
    epoch_list2 = [1, 2, 3, 4, 5, 6, 7]
    expected_output2 = {
        'copy_number': 6,
        'epoch_created': 1,
        'child': {
            'copy_number': 4,
            'epoch_created': 2,
            'child': {
                'copy_number': 3,
                'epoch_created': 3
            },
            'complement': {
                'copy_number': 1,
                'epoch_created': 4
            }
        },
        'complement': {
            'copy_number': 2,
            'epoch_created': 5,
            'child': {
                'copy_number': 1,
                'epoch_created': 6
            },
            'complement': {
                'copy_number': 1,
                'epoch_created': 7
            }
        }
    }
    assert add_epoch_created(tree2, epoch_list2) == expected_output2, f"Error: {add_epoch_created(tree2, epoch_list2)}"

    # Test 3: Tree with fewer nodes than epochs
    tree3 = {'copy_number': 2, 'child': {'copy_number': 1}, 'complement': {'copy_number': 1}}
    epoch_list3 = [1, 2]
    try:
        add_epoch_created(tree3, epoch_list3)
    except ValueError as e:
        assert str(e) == "Epoch list is shorter than the number of nodes in the tree. Tree: {'copy_number': 2, 'child': {'copy_number': 1, 'epoch_created': 2}, 'complement': {'copy_number': 1}, 'epoch_created': 1}, Epoch List: [1, 2]", f"Error: {e}"
    else:
        raise AssertionError("Expected ValueError but none was raised.")

    # Test 4: Tree with fewer epochs than nodes
    tree4 = {'copy_number': 2, 'child': {'copy_number': 1}, 'complement': {'copy_number': 1}}
    epoch_list4 = [1, 2, 3, 4]
    try:
        add_epoch_created(tree4, epoch_list4)
    except ValueError as e:
        assert str(e) == "Epoch list is longer than the number of nodes in the tree. Tree: {'copy_number': 2, 'child': {'copy_number': 1, 'epoch_created': 2}, 'complement': {'copy_number': 1, 'epoch_created': 3}, 'epoch_created': 1}, Epoch List: [1, 2, 3, 4]", f"Error: {e}"
    else:
        raise AssertionError("Expected ValueError but none was raised.")

test_add_epoch_created()


def update_to_epoch_created(tree, parent_epoch=None):
    """
    Recursively update the 'epoch_created' field of a tree and its children.
    
    :param tree: A dictionary representing a tree with 'child' and 'complement' keys.
    :param parent_epoch: The epoch value to update the tree with. If None, the tree's own value is used.
    :return: The updated tree.
    """
    def _update_subtree(subtree_key):
        if subtree_key in tree and tree[subtree_key] is not None:
            update_to_epoch_created(tree[subtree_key], tree['epoch_created'])

    _update_subtree('child')
    _update_subtree('complement')

    if parent_epoch is not None:
        tree['epoch_created'] = parent_epoch

    return tree


def test_update_to_epoch_created():
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

    updated_tree = update_to_epoch_created(tree)
    expected_tree = {
        'epoch_created': 1,
        'child': {
            'epoch_created': 1,
            'child': None,
            'complement': None
        },
        'complement': {
            'epoch_created': 1,
            'child': None,
            'complement': None
        }
    }
    
    try:
        assert updated_tree == expected_tree
    except AssertionError:
        pretty_print(f"Assertion failed: expected {expected_tree}, got {updated_tree}")
        raise


test_update_to_epoch_created()

def convert_to_dict_tree(tree):
    """
    Convert a tree data structure in the form [value, [child1], [child2]] or (value, (child1), (child2))
    into a dictionary-like data structure in the form
    {'copy_number': value, 'child': child1, 'complement': child2}.
    
    :param tree: list or tuple representing the tree data structure
    :return: dictionary-like data structure
    """
    if isinstance(tree, (list, tuple)) and len(tree) == 3:
        copy_number, child1, child2 = tree

        # Check if child1 and child2 are trees (list or tuple of length 3 or 1)
        # If they are, then recursively process them
        if isinstance(child1, (list, tuple)) and len(child1) in [1, 3]:
            child1 = convert_to_dict_tree(child1)
        else:
            raise ValueError(f"Child node is not in the expected format. Got: {child1}")
        
        if isinstance(child2, (list, tuple)) and len(child2) in [1, 3]:
            child2 = convert_to_dict_tree(child2)
        else:
            raise ValueError(f"Complement node is not in the expected format. Got: {child2}")

        return {'copy_number': copy_number, 'child': child1, 'complement': child2}

    elif isinstance(tree, (list, tuple)) and len(tree) == 1:
        return {'copy_number': tree[0]}

    else:
        raise ValueError(f"Input tree is not in the expected format. Got: {tree}")


def convert_to_dict_tree(tree):
    """
    Convert a tree data structure in the form [value, [child1], [child2]] or (value, (child1), (child2))
    into a dictionary-like data structure in the form
    {'copy_number': value, 'child': child1, 'complement': child2}.
    
    :param tree: list or tuple representing the tree data structure
    :return: dictionary-like data structure
    """
    
    # Check for single element in list or tuple
    if isinstance(tree, (list, tuple)) and len(tree) == 1:
        result = {'copy_number': tree[0]}
    
    # Check for three elements in list or tuple
    elif isinstance(tree, (list, tuple)) and len(tree) == 3:
        copy_number, child1, child2 = tree
        result = {'copy_number': copy_number, 
                'child': convert_to_dict_tree(child1), 
                'complement': convert_to_dict_tree(child2)}
    
    # Raise error for invalid input
    else:
        raise ValueError(f"Input tree is not in the expected format. Got: {tree}")


    return result


def test_convert_to_dict_tree():
    # Test 1: Tree with one node
    tree1 = [3]
    expected_output1 = {'copy_number': 3}
    assert convert_to_dict_tree(tree1) == expected_output1, f"Error: {convert_to_dict_tree(tree1)}"

    # Test 2: Tree with two levels
    tree2 = [4, [2, [1]], [2, [1]]]
    expected_output2 = {'copy_number': 4, 'child': {'copy_number': 2, 'child': {'copy_number': 1}}, 'complement': {'copy_number': 2, 'child': {'copy_number': 1}}}
    assert convert_to_dict_tree(tree2) == expected_output2, f"Error: {convert_to_dict_tree(tree2)}"

    # Test 3: More complex tree
    tree3 = [6, [4, [3], [1]], [2, [1], [1]]]
    expected_output3 = {'copy_number': 6, 'child': {'copy_number': 4, 'child': {'copy_number': 3}, 'complement': {'copy_number': 1}}, 'complement': {'copy_number': 2, 'child': {'copy_number': 1}, 'complement': {'copy_number': 1}}}
    assert convert_to_dict_tree(tree3) == expected_output3, f"Error: {convert_to_dict_tree(tree3)}"

    # Test 4: Tree with invalid input
    tree4 = [1, [2, 3], 4]
    try:
        convert_to_dict_tree(tree4)
    except ValueError as e:
        assert str(e) == "Input tree is not in the expected format.", f"Error: {e}"
    else:
        raise AssertionError("Expected ValueError but none was raised.")


def order_tree_keys_alphabetically(tree):
    """
    Recursively orders the keys of a given dictionary (tree) alphabetically.
    The function processes nested dictionaries with the keys 'child' and 'complement'.
    
    Args:
    tree (dict): The dictionary (tree) with keys to be ordered alphabetically.
    
    Returns:
    dict: A new dictionary with keys ordered alphabetically.
    """
    ordered_tree = {}
    
    for key in sorted(tree.keys(), reverse=True):
        if key in ['child', 'complement']:
            ordered_tree[key] = order_tree_keys_alphabetically(tree[key]) if tree[key] else None
        else:
            ordered_tree[key] = tree[key]
    
    return ordered_tree


def test_order_tree_keys_alphabetically():
    test_tree = {
        "B": "value B",
        "A": "value A",
        "complement": {
            "D": "value D",
            "C": "value C",
            "child": {
                "F": "value F",
                "E": "value E"
            }
        }
    }

    expected_result = {
        "A": "value A",
        "B": "value B",
        "complement": {
            "C": "value C",
            "D": "value D",
            "child": {
                "E": "value E",
                "F": "value F"
            }
        }
    }

    result = order_tree_keys_alphabetically(test_tree)
    assert result == expected_result, f"Expected {expected_result}, but got {result}"

    empty_tree = {}
    expected_empty_tree = {}
    result = order_tree_keys_alphabetically(empty_tree)
    assert result == expected_empty_tree, f"Expected {expected_empty_tree}, but got {result}"



test_order_tree_keys_alphabetically()



def convert_CN_tree_and_epoch_list_to_dict_tree(CN_tree, epoch_list):
    """
    Convert the given CN_tree and epoch_list into a dictionary tree.
    
    :param CN_tree: A nested data structure representing a CN tree.
    :param epoch_list: A list of epoch values.
    :return: A dictionary tree created from the CN_tree and epoch_list.
    """
    epoch_list = list(epoch_list)
    dict_tree = convert_to_dict_tree(CN_tree)
    dict_tree = add_epoch_created(dict_tree, epoch_list)
    return dict_tree


def filter_tree(tree, keys_to_keep):
    """
    Filter the given tree, keeping only the specified keys.
    
    :param tree: A dictionary representing a tree structure.
    :param keys_to_keep: A list of keys to keep in the filtered tree.
    :return: A new tree with only the specified keys and their values.
    """
    filtered_tree = {key: tree[key] for key in keys_to_keep if key in tree}
    
    if 'child' in tree and tree['child'] is not None:
        filtered_tree['child'] = filter_tree(tree['child'], keys_to_keep)
        
    if 'complement' in tree and tree['complement'] is not None:
        filtered_tree['complement'] = filter_tree(tree['complement'], keys_to_keep)
        
    return filtered_tree


def test_convert_CN_tree_and_epoch_list_to_dict_tree():
    CN_tree = ("A", ("B", None, None), ("C", None, None))
    epoch_list = [1, 2, 3]
    expected_output = {
        'copy_number': 'A',
        'child': {
            'copy_number': 'B',
            'child': None,
            'complement': None,
            'epoch_created': 2
        },
        'complement': {
            'copy_number': 'C',
            'child': None,
            'complement': None,
            'epoch_created': 3
        },
        'epoch_created': 1
    }

    result = convert_CN_tree_and_epoch_list_to_dict_tree(CN_tree, epoch_list)
    assert result == expected_output, f"Expected {expected_output}, but got {result}"


def test_filter_tree():
    tree = {
        "name": "A",
        "epoch_created": 1,
        "child": {
            "name": "B",
            "epoch_created": 2,
            "child": None,
            "complement": None,
        },
        "complement": {
            "name": "C",
            "epoch_created": 3,
            "child": None,
            "complement": None,
        },
    }
    keys_to_keep = ["name", "complement"]
    expected_output = {
        'name': 'A',
        'complement': {
            'name': 'C',
            'complement': None
        },
        'child': {
            'name': 'B',
            'complement': None
        }
    }

    result = filter_tree(tree, keys_to_keep)
    assert result == expected_output, f"Expected {expected_output}, but got {result}"


test_filter_tree()



def create_epoch_index(tree, key_from, key_to):
    """
    The function replaces the key called 'key_from' in the dictionary 'tree' with a new key called 'key_to',
    keeping the original value.

    :param tree: A dictionary representing a tree node.
    :param key_from: The original key to be replaced.
    :param key_to: The new key that will replace the original key.
    :return: The modified tree with the key replaced.
    """
    if key_from in tree:
        tree[key_to] = tree[key_from]
        del tree[key_from]

    for child_key in ['child', 'complement']:
        if child_key in tree and tree[child_key] is not None:
            tree[child_key] = create_epoch_index(tree[child_key], key_from, key_to)

    return tree


def test_create_epoch_index():
    tree = {
        'a': 1,
        'child': {'b': 2, 'child': {'c': 3}},
        'complement': {'d': 4}
    }
    expected_result = {
        'x': 1,
        'child': {'b': 2, 'child': {'c': 3}},
        'complement': {'d': 4}
    }
    result = create_epoch_index(tree, 'a', 'x')

    assert result == expected_result, f"Expected: {expected_result}, got: {result}"

    tree = {
        'a': 1,
        'child': {'b': 2, 'child': {'c': 3}},
        'complement': {'d': 4}
    }
    expected_result = {
        'a': 1,
        'child': {'b': 2, 'child': {'z': 3}},
        'complement': {'d': 4}
    }
    result = create_epoch_index(tree, 'c', 'z')

    assert result == expected_result, f"Expected: {expected_result}, got: {result}"


test_create_epoch_index()


def pretty_print_tree(tree):
    """
    Pretty-prints a nested dictionary representing a tree.

    Args:
        tree (dict): The tree dictionary to be pretty-printed.
    """
    print(json.dumps(tree, indent=4, sort_keys=True))



def print_summary(total_nodes, num_chrom_with_correct_CN, num_chrom_with_correct_CN_and_epoch_created, average_distance_from_truth_of_epoch_created):
    print(f"Total nodes: {total_nodes}")
    print(f"Number of chromosomes with correct copy numbers: {num_chrom_with_correct_CN}")
    print(f"Number of chromosomes with correct copy numbers and epochs created: {num_chrom_with_correct_CN_and_epoch_created}")
    print(f"Average distance from truth of epoch created: {average_distance_from_truth_of_epoch_created}")



def load_simulation_data(test_case, simulation_filename):
    """
    This function opens a file containing simulation data and returns the simulation results.

    Args:
    test_case (str): A string specifying which test case to load.
    simulation_filename (str): A string specifying the name of the simulation.

    Returns:
    dict: A dictionary containing the results of the simulation.
    """
    with open(f'{SIMULATIONS_FILE_FOLDER}/{simulation_filename}_{test_case}.pickle', 'rb') as f:
        SS = pickle.load(f)
    return SS


def save_simulation_data(test_case, simulation_filename, data):
    """
    This function saves simulation data to a file.

    Args:
    test_case (str): A string specifying which test case to save.
    simulation_filename (str): A string specifying the name of the simulation.
    data (dict): A dictionary containing the results of the simulation.
    """
    filename = f'{SIMULATIONS_FILE_FOLDER}/{simulation_filename}_{test_case}.pickle'
    filename_smaller = f'{SIMULATIONS_FILE_FOLDER}/{simulation_filename}_{test_case}_smaller.pickle'

    if os.path.exists(filename):
        with open(filename, 'rb') as f:
            old_data = pickle.load(f)
        old_data.update(data)
        data = old_data

    # Save the full version
    with open(filename, 'wb') as f:
        pickle.dump(data, f)

    # Create a smaller version by removing certain keys
    smaller_data = data.copy()
    if 'likelihoods' in smaller_data:
        del smaller_data['likelihoods']
    if 'simulated_chromosomes' in smaller_data:
        del smaller_data['simulated_chromosomes']

    # Save the smaller version
    with open(filename_smaller, 'wb') as f:
        pickle.dump(smaller_data, f)


def analyze_sorted_results(SS):
    print("Do some analysis")
    for solution in SS["solutions"]:
        total_truth_nodes = 0
        total_solution_nodes = 0
        num_chrom_with_correct_CN = 0
        num_chrom_with_correct_CN_and_epoch_created = 0
        average_distance_from_truth_of_epoch_created = 0
        SS["total_nodes"] = {}

        for chrom in range(23):
            # are trees same by epoch and copy number


            # count the total nodes:
            solution[chrom]["total_nodes"] = count_nodes(solution[chrom]["dict_tree"])
            total_solution_nodes += solution[chrom]["total_nodes"]

            SS["total_nodes"][chrom] = count_nodes(SS["simplified_truth_trees"][chrom])
            total_truth_nodes += SS["total_nodes"][chrom]

            # count the number with the correct CN:
            solution[chrom]["correct_CN"] = analyze_tree_copy_number(SS["simplified_truth_trees"][chrom], solution[chrom]["dict_tree"])
            num_chrom_with_correct_CN += solution[chrom]["correct_CN"]

            # count the number with the correct CN and epoch created:
            solution[chrom]["correct_CN_and_epoch_created"] = analyze_tree_properties(SS["simplified_truth_trees"][chrom], solution[chrom]["dict_tree"])
            num_chrom_with_correct_CN_and_epoch_created += solution[chrom]["correct_CN_and_epoch_created"]

            # find the average epoch distance for tree of correct CN (incorrect CN trees get infinite distance)
            solution[chrom]["error_in_epoch_created_estimate"] = calculate_tree_distance(SS["simplified_truth_trees"][chrom], solution[chrom]["dict_tree"], total_nodes)
            average_distance_from_truth_of_epoch_created += solution[chrom]["error_in_epoch_created_estimate"]


    print_summary(total_nodes, num_chrom_with_correct_CN,
                  num_chrom_with_correct_CN_and_epoch_created, average_distance_from_truth_of_epoch_created)

    for result_index, res in enumerate(sorted_results):
        new_result = compile_new_result(res, SS, pre_est, mid_est, post_est)

        process_new_result(new_result, SS)

        file_name = f'{SIMULATIONS_FILE_FOLDER}/{simulation_filename}_{test_case}.pickle'
        all_results = load_simulation_data(test_case, simulation_filename)
        all_results[str(result_index)] = new_result
        write_simulation_data(test_case, simulation_filename, all_results)

def get_all_trees_ready_for_comparison(SS):
    #print("get trees ready for comparison")

    get_truth_trees_ready_for_comparison(SS)
    #print("truth trees are ready")

    get_estimated_trees_ready_for_comparison(SS)
    #print("estimated trees are ready")


def get_truth_trees_ready_for_comparison(SS):
    #pretty_print("simplified simulated tree")

    if "simplified_truth_trees" not in SS:
        SS["simplified_truth_trees"] = {}

    assert(sorted(list(SS["truth_trees"].keys())) == list(range(23))) # assert it has all the keys available

    for chrom in SS["truth_trees"]:
        #print("Current chromosome: ", chrom)
        SS["simplified_truth_trees"][chrom] = filter_tree(SS["truth_trees"][chrom],keys_to_keep = ["copy_number","epoch_created"])
        #print("Newly added tree to simplified_truth_trees: ", SS["simplified_truth_trees"][chrom])

        SS["simplified_truth_trees"][chrom] = create_epoch_index(SS["simplified_truth_trees"][chrom], key_from="epoch_created", key_to="epoch_index")
        #print("Trees after adding epoch index: ",  SS["simplified_truth_trees"][chrom])

        SS["simplified_truth_trees"][chrom] = order_tree_keys_alphabetically(SS["simplified_truth_trees"][chrom])
        #print("Tree after ordering keys alphabetically: ", SS["simplified_truth_trees"][chrom])

        SS["simplified_truth_trees"][chrom] = sort_tree_by_copynumber(SS["simplified_truth_trees"][chrom])
        #print("Tree after ordering keys by copynumber: ", SS["simplified_truth_trees"][chrom])
        #print("\n"*5)


def get_estimated_trees_ready_for_comparison(SS):
    for solution in SS["solutions"]:
        for chrom in range(23):
            assert(chrom in solution)
            solution[chrom]["dict_tree"] = convert_CN_tree_and_epoch_list_to_dict_tree(solution[chrom]["tree"], solution[chrom]["epochs_created"])

            solution[chrom]["dict_tree"] = create_epoch_index(solution[chrom]["dict_tree"], key_from="epoch_created", key_to="epoch_index")

            solution[chrom]["dict_tree"] = order_tree_keys_alphabetically(solution[chrom]["dict_tree"])

            solution[chrom]["dict_tree"] = sort_tree_by_copynumber(solution[chrom]["dict_tree"])

            solution[chrom]["dict_tree"]['epoch_index'] = 0


def sort_simulation_results_by_likelihood(solutions):
    """
    This function sorts simulation results based on 'best_loglik'.
    
    Args:
    results (list): A list of dictionaries containing simulation results.

    Returns:
    list: A sorted list of dictionaries.
    """
    return sorted(solutions, key=lambda x: x['est_neg_loglik'], reverse=True)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Run all processing scripts with the given test case and optional parameters."
    )

    parser.add_argument(
        "-t", 
        "--test_case", 
        type=int, 
        help="The name or number of the test case you want to process."
    )

    parser.add_argument(
        "-f", 
        "--simulation_filename", 
        required=True, 
        help="Name of the simulation."
    )

    return parser.parse_args()


def print_tree(tree):
    """
    Takes a dictionary tree and prints a neatly formatted list representation.
    
    :param tree: Dictionary tree to be printed.
    """
    converted_tree = convert_dict_tree_to_list(tree)
    
    def print_tree_recursively(tree, indent=0):
        """
        Helper function to recursively print a tree.
        
        :param tree: List representation of the tree.
        :param indent: Current indentation level (default: 0).
        """
        print('\t' * indent, tree[0])
        for subtree in tree[1:]:
            if subtree is not None:
                print_tree_recursively(subtree, indent + 1)

    print_tree_recursively(converted_tree)


def print_solution_results(SS):
    """
    Prints the solution results for intuition. For each chromosome in the solution,
    this function prints the dictionary representation of the solution tree and the corresponding
    truth tree from the simplified truth trees.

    :param SS: A dictionary containing the solutions and the simplified truth trees.
    """
    print("Print results for intuition:")
    for idx, solution in enumerate(SS["solutions"]):
        print(f"\nSolution {idx}:")

        # Print non-integer keys:
        for key, value in solution.items():
            if not isinstance(key, int):
                print(f"{key}: {value}")

        for chrom in range(23):
            print("\n\tSolution Tree for Chromosome", chrom)
            print(solution[chrom]["dict_tree"])
            print_tree(solution[chrom]["dict_tree"])
            print("\n\tTruth Tree for Chromosome", chrom)
            print(SS["simplified_truth_trees"][chrom])
            print_tree(SS["simplified_truth_trees"][chrom])



def print_similarity_results(SS: Dict):
    """
    Print the similarity results computed in each solution along with comparison results of relative timings.

    :param SS: A dictionary containing the solutions, the computed similarity results, and comparison results.
    """
    print("Print similarity results:")

    # Set the display options.
    pd.set_option('display.max_columns', None) 
    pd.set_option('display.max_rows', None) 
    pd.set_option('display.max_colwidth', None)
    pd.set_option('display.width', 200)  # default is 80


    for idx, solution in enumerate(SS["solutions"]):
        print(f"\nSolution {idx}:")

        # Initialize an empty DataFrame for each solution
        df = pd.DataFrame(columns=['Copy Number Distance Function', 'Epoch Created Distance Function', 'Normalising Constant', 'Similarity Score'])

        # Print the keys related to similarity results
        for key, value in solution.items():
            if isinstance(key, str):
                if any(distance_function in key for distance_function in ["ZeroDistance", "AbsoluteDifference", "SquaredDistance", "InfiniteDistance"]):
                    distance_function_1, distance_function_2, constant = key.split('_')

                    # Append a new row to the DataFrame
                    new_row = pd.DataFrame([{
                        'Copy Number Distance Function': distance_function_1,
                        'Epoch Created Distance Function': distance_function_2,
                        'Normalising Constant': constant,
                        'Similarity Score': value
                    }])
                    df = pd.concat([df, new_row], ignore_index=True)

        solution["distance_dataframe"] = df

        # Print the DataFrame
        df = df.sort_values(by=['Normalising Constant', 'Copy Number Distance Function'])
        print(df)

        solution["correct_path"] = SS['pre'] == solution[0]['pre'] and SS['mid'] == solution[0]['mid'] and SS['post'] == solution[0]['post']
         
        solution_data = {
            'best_neg_loglik': [solution['est_neg_loglik']],
            'best_p_up': [solution['est_p_up']],
            'best_p_down': [solution['est_p_down']],
            'best_plambda': [solution['est_plambda']],
            'pre': [solution[0]['pre']],
            'mid': [solution[0]['mid']],
            'post': [solution[0]['post']],
            'total_time': [solution['execution_times']['total_time']],
            'time_get_all_trees_and_timings': [solution['execution_times']['get_all_trees_and_timings']],
            'time_timing_struct_to_all_structures': [solution['execution_times']['timing_struct_to_all_structures']],
            'time_find_BP_and_SNV_loglik': [solution['execution_times']['timing_struct_to_all_structures']],
            'counts_<': [solution['counts_<']],
            'counts_<=': [solution['counts_<=']],
            'correct_path': solution["correct_path"]
        }

        df = pd.DataFrame(solution_data)
        print(df)
        solution["metadata_dataframe"] = df


############
############
############
############
############
############
class DistanceFunction(ABC):
    """
    Abstract base class for distance functions. 
    Defines an interface for a function that takes an integer difference 
    and returns a float.
    """
    
    @abstractmethod
    def __call__(self, difference: int) -> float:
        """
        Abstract method that will be implemented in the derived classes.

        :param difference: An integer to be used in the distance function.
        :return: A float as a result of the distance function.
        """
        pass


class ZeroDistance(DistanceFunction):
    """
    Zero distance function class. Always returns zero regardless of the input.
    """
    
    def __call__(self, difference: int) -> float:
        assert isinstance(difference, int), f"Invalid input type. Expected int but got {type(difference)}"
        return 0.0

class AbsoluteDifference(DistanceFunction):
    """
    Euclidean distance function class.
    """
    def __call__(self, difference: int) -> float:
        assert isinstance(difference, int), f"Invalid input type. Expected int but got {type(difference)}"
        return abs(float(difference))


class SquaredDistance(DistanceFunction):
    """
    Squared distance function class.
    """
    def __call__(self, difference: int) -> float:
        assert isinstance(difference, int), f"Invalid input type. Expected int but got {type(difference)}"
        return float(difference ** 2)


class InfiniteDistance(DistanceFunction):
    """
    Infinite distance function class.
    """
    def __call__(self, difference: int) -> float:
        assert isinstance(difference, int), f"Invalid input type. Expected int but got {type(difference)}"
        if difference != 0:
            return float('inf')
        else:
            return 0


def test_distance_function(distance_function: DistanceFunction, difference: int, expected_output: float):
    """
    Test function for the DistanceFunction class and its subclasses.

    :param distance_function: An instance of a DistanceFunction subclass.
    :param difference: An integer to be used in the distance function.
    :param expected_output: The expected output of the distance function.
    """
    assert isinstance(difference, int), f"Invalid input type. Expected int but got {type(difference)}"
    assert isinstance(expected_output, float), f"Invalid expected output type. Expected float but got {type(expected_output)}"
    
    result = distance_function(difference)
    assert result == pytest.approx(expected_output), f"Expected {expected_output}, but got {result}"

def run_tests():
    # Testing ZeroDistance class
    zero_distance = ZeroDistance()
    test_distance_function(zero_distance, 5, 0.0)
    test_distance_function(zero_distance, 0, 0.0)

    # Testing AbsoluteDifference class
    euclidean_distance = AbsoluteDifference()
    test_distance_function(euclidean_distance, 5, 5.0)
    test_distance_function(euclidean_distance, 0, 0.0)

    # Testing SquaredDistance class
    squared_distance = SquaredDistance()
    test_distance_function(squared_distance, 2, 4.0)
    test_distance_function(squared_distance, 0, 0.0)



run_tests()



def calculate_difference(tree1, tree2, distance_function):
    return distance_function(abs(tree1 - tree2))

def compute_similarity(
    tree1: Union[Dict, None], 
    tree2: Union[Dict, None], 
    copy_number_distance: DistanceFunction, 
    epoch_index_distance: DistanceFunction, 
    similarity_divisor: int
) -> float:
    # Handling None cases for tree1 and tree2
    if tree1 is None and tree2 is None:
        return 0.0

    elif tree1 is None or tree2 is None:
        if tree1 is None:
            tree1 = {'copy_number': 0, 'epoch_index': 0}
        else:
            assert(tree2 is None)
            tree2 = {'copy_number': 0, 'epoch_index': 0}

        return calculate_difference(
            tree1['copy_number'], 
            tree2['copy_number'], 
            copy_number_distance
        )
    
    else:
        epoch_index_difference = calculate_difference(
            tree1['epoch_index'], 
            tree2['epoch_index'], 
            epoch_index_distance
        )

        left_similarity = compute_similarity(
            tree1.get('complement', None), 
            tree2.get('complement', None), 
            copy_number_distance, 
            epoch_index_distance, 
            similarity_divisor
        )

        right_similarity = compute_similarity(
            tree1.get('child', None), 
            tree2.get('child', None), 
            copy_number_distance, 
            epoch_index_distance, 
            similarity_divisor
        )

        #print("*****",tree1,tree2,epoch_index_difference, left_similarity, right_similarity)

        total_similarity = (
            epoch_index_difference 
            + left_similarity 
            + right_similarity
        ) / similarity_divisor

        return total_similarity



def test_calculate_difference():
    calc = NodeSimilarityCalculator()
    tree1 = {"copy_number": 5, "epoch_index": 10}
    tree2 = {"copy_number": 3, "epoch_index": 8}

    copy_number_distance = lambda x: x
    epoch_index_distance = lambda x: x
    assert calc._calculate_difference(tree1, tree2, copy_number_distance, epoch_index_distance) == (2, 2)

def test_compute_similarity():
    calc = NodeSimilarityCalculator()
    tree1 = {"copy_number": 5, "epoch_index": 10}
    tree2 = {"copy_number": 3, "epoch_index": 8}

    copy_number_distance = lambda x: x
    epoch_index_distance = lambda x: x
    assert calc.compute_similarity(tree1, tree2, copy_number_distance, epoch_index_distance) == (0.5, 1.0)


def compute_solution_similarity(
    solutions: Dict,
    node_distance_functions: List[Type[DistanceFunction]]
) -> None:
    """
    Compute similarity for solutions with various distance and similarity functions.

    Args:
        solutions: Solutions dictionary containing trees to compare.
        node_distance_functions: List of distance function classes for node comparisons.

    Returns: None. The function modifies the solutions dictionary in-place.
    """
    # Create distance function pairs
    pairs = [(x, y) for y in node_distance_functions for x in node_distance_functions]

    # Iterate over solutions
    for solution in solutions["solutions"]:
        for node_distance_function_pair in pairs:
            for normalising_constant in [1,4]:

                # Create instances of distance function classes
                copy_number_distance = node_distance_function_pair[0]()
                epoch_created_distance = node_distance_function_pair[1]()

                # Initialize total similarity
                total_similarity = 0

                for chrom in range(23):
                    # Compute similarity
                    similarity = compute_similarity(
                        solution[chrom]["dict_tree"],
                        solutions["simplified_truth_trees"][chrom],
                        copy_number_distance,
                        epoch_created_distance,
                        normalising_constant
                    )

                    # Add to total similarity
                    total_similarity += similarity

                # Save result
                key = f"{copy_number_distance.__class__.__name__}_{epoch_created_distance.__class__.__name__}_{normalising_constant}"
                solution[key] = total_similarity



def extract_timing_values(node: Dict, epochs: List[Tuple[int]], is_root: bool = True) -> None:
    """
    Recursive function to extract timing values from a tree structure.

    :param node: The current node in the tree.
    :param epochs: The list of timing values collected so far.
    :param is_root: A flag to indicate if the node is the root node.
    """
    if not is_root:  # Only append if it's not the root
        epochs.append((node['epoch_index'],))
        
    if 'complement' in node:
        extract_timing_values(node['complement'], epochs, is_root=False)
    if 'child' in node:
        extract_timing_values(node['child'], epochs, is_root=False)


def align_timing_values(truth_epochs: List[Tuple[int]], solution_epochs: List[Tuple[int]]) -> Tuple[List[Tuple[int]], List[Tuple[int]]]:
    """
    Aligns the timing value lists by padding the shorter one with the last element of the longer list.

    :param truth_epochs: List of truth timing values.
    :param solution_epochs: List of solution timing values.
    :return: The aligned lists of truth and solution timing values.
    """
    if len(truth_epochs) < len(solution_epochs):
        truth_epochs += [truth_epochs[-1]] * (len(solution_epochs) - len(truth_epochs))
    elif len(solution_epochs) < len(truth_epochs):
        solution_epochs += [solution_epochs[-1]] * (len(truth_epochs) - len(solution_epochs))
    return truth_epochs, solution_epochs

def compare_relative_timing(SS, solution: Dict, operator: str, is_root: bool = True) -> Tuple[int, int, int]:
    """
    Compares relative timing of nodes across multiple chromosomes using a specified operator.

    :param truth: The truth tree structure.
    :param solution: The solution tree structure.
    :param operator: The operator to use for comparison ("<" or "<=").
    :return: The count of correct, incorrect and missing comparisons.
    """
    count_true = 0
    count_false = 0
    count_missing = 0

    for chromosome in SS["simplified_truth_trees"]:
        truth_epochs = []
        solution_epochs = []

        extract_timing_values(SS["simplified_truth_trees"][chromosome], truth_epochs, is_root)
        extract_timing_values(solution[chromosome]["dict_tree"], solution_epochs, is_root)

        truth_epochs, solution_epochs = align_timing_values(truth_epochs, solution_epochs)

        for i in range(len(truth_epochs)):
            for j in range(i + 1, len(truth_epochs)):
                if operator == "<":
                    if truth_epochs[i] < truth_epochs[j] and solution_epochs[i] < solution_epochs[j]:
                        count_true += 1
                    elif truth_epochs[i] < truth_epochs[j] and solution_epochs[i] >= solution_epochs[j]:
                        count_false += 1
                elif operator == "<=":
                    if truth_epochs[i] <= truth_epochs[j] and solution_epochs[i] <= solution_epochs[j]:
                        count_true += 1
                    elif truth_epochs[i] <= truth_epochs[j] and solution_epochs[i] > solution_epochs[j]:
                        count_false += 1

                if solution_epochs[i] is None or solution_epochs[j] is None:
                    count_missing += 1

    return count_true, count_false, count_missing


def add_relative_timing_comparison(SS: Dict, is_root: bool = True) -> Dict:
    """
    Applies the compare_relative_timing function to each solution in SS["solutions"] and stores the results in new keys.

    :param SS: The original structure of solutions and truth.
    :return: The updated structure with added comparison results.
    """
    for solution in SS["solutions"]:
        for operator in ["<", "<="]:
            count_true, count_false, count_missing = compare_relative_timing(SS, solution, operator, is_root)
            solution["counts_" + operator] = {
                "count_true": count_true,
                "count_false": count_false,
                "count_missing": count_missing
            }
    return SS



############
############
############
############
############
############
def add_ev_strings_and_counts_to_dicts(SS):
    truth_str = get_ev_string(SS['pre'], SS['mid'], SS['post'])
    true_GD_count = count_genome_doublings(truth_str)
    for d in SS["solutions"]:
        d['ev_str'] = truth_str
        d['est_pre'] = d[0]['pre']
        d['est_mid'] = d[0]['mid']
        d['est_post'] = d[0]['post']
        d['pre'] = SS['pre']
        d['mid'] = SS['mid']
        d['post'] = SS['post']
        d['p_up'] = SS['p_up']
        d['p_down'] = SS['p_down']
        d['plambda'] = SS['rate']

        d['est_ev_str'] = get_ev_string(d['est_pre'], d['est_mid'], d['est_post'])
        d['genome_doublings'] = true_GD_count
        d['est_genome_doublings'] = count_genome_doublings(d['est_ev_str'])


def add_aic_to_dicts(SS):
    for d in SS["solutions"]:
        num_parameters = d['est_genome_doublings'] + 3
        neg_log_likelihood = d['est_neg_loglik']
        aic = 2 * num_parameters + 2 * neg_log_likelihood
        d['AIC'] = aic


def sort_dicts_by_worst_aic(SS):
    SS["solutions"] = sorted(SS["solutions"], key=lambda x: -x['AIC'])

def count_genome_doublings(ev_string):
    return ev_string.count("G")

def process_further(SS):
    add_ev_strings_and_counts_to_dicts(SS)
    add_aic_to_dicts(SS)
    sort_dicts_by_worst_aic(SS)


def get_size(obj):
    size = sys.getsizeof(obj)
    if isinstance(obj, dict):
        size += sum(get_size(value) for value in obj.values())
    elif isinstance(obj, list):
        size += sum(get_size(item) for item in obj)
    return size

def print_dicts_memory(ss, indent=0, level=1, max_level=4):
    if level > max_level:
        return

    for key, value in ss.items():
        size = get_size(value)
        if isinstance(value, dict):
            print(' ' * indent + str(key) + ': dict of total size ' + str(size) + ' bytes')
            print_dicts_memory(value, indent + 4, level + 1, max_level)
        elif isinstance(value, list):
            print(' ' * indent + str(key) + ': list of total size ' + str(size) + ' bytes')
        else:
            print(' ' * indent + str(key) + ': ' + str(size) + ' bytes')


def main(test_case: str, simulation_filename: str) -> None:
    """
    The main function that processes a simulation test case and analyzes estimated trees.

    Args:
        test_case (str): A string specifying which test case to process.
        simulation_filename (str): The filename containing the simulation data.

    """
    SS = load_simulation_data(test_case, simulation_filename)
    SS["solutions"] = sort_simulation_results_by_likelihood(SS["solutions"])
    get_all_trees_ready_for_comparison(SS)
    print_solution_results(SS)

    # List of node distance functions
    node_distance_functions = [ZeroDistance, AbsoluteDifference, SquaredDistance, InfiniteDistance]

    # Call the function
    compute_solution_similarity(SS, node_distance_functions)
    add_relative_timing_comparison(SS)
    print_similarity_results(SS)

    process_further(SS)

    save_simulation_data(test_case, simulation_filename, SS)

    
    print_dicts_memory(SS)

if __name__ == "__main__":
    args = parse_arguments()
    test_case = args.test_case
    simulation_filename = args.simulation_filename
    main(test_case, simulation_filename)


