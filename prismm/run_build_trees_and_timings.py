import pickle as pkl
import scipy
import numpy as np
import copy
import re
from typing import Dict
from constants import *
from general_functions import *
from run_simulation import generate_dirichlet_probability
import scipy.optimize as opt
#from scipy.optimize import minimize_scalar
import sys
import argparse
import logging
import pickle
import os



logging.basicConfig(level = logging.INFO)
###### STEP 3; calculate the SNV multiplicities of each chromosome
######
######
######
######

def simulated_chromosomes_to_SNV_counts(simulated_chromosomes):
    """
    Count the copy number of each SNV in the simulated chromosomes.

    :param simulated_chromosomes: A dictionary containing the simulated chromosomes.
    :return: A dictionary with the count of copy numbers for each SNV in the simulated chromosomes.
    """

    SNV_copy_counter = {}

    for chrom_type in simulated_chromosomes:
        if chrom_type not in SNV_copy_counter:
            SNV_copy_counter[chrom_type] = {}

        for chrom in simulated_chromosomes[chrom_type]:
            if chrom["dead"]:
                continue

            for SNV in chrom["SNVs"]:
                UI = SNV["unique_identifier"]

                if UI not in SNV_copy_counter[chrom_type]:
                    SNV_copy_counter[chrom_type][UI] = 1
                else:
                    SNV_copy_counter[chrom_type][UI] += 1

    return SNV_copy_counter


# Test case for simulated_chromosomes_to_SNV_counts
def test_simulated_chromosomes_to_SNV_counts():
    simulated_chromosomes = {
        "A": [{"dead": False, "SNVs": [{"unique_identifier": "1"}, {"unique_identifier": "2"}]},
              {"dead": True, "SNVs": [{"unique_identifier": "1"}, {"unique_identifier": "3"}]}],
        "B": [{"dead": False, "SNVs": [{"unique_identifier": "2"}, {"unique_identifier": "3"}]},
              {"dead": False, "SNVs": [{"unique_identifier": "1"}, {"unique_identifier": "4"}]}],
        "C": [{"dead": False, "SNVs": [{"unique_identifier": "2"}, {"unique_identifier": "3"}]},
              {"dead": False, "SNVs": [{"unique_identifier": "1"}, {"unique_identifier": "2"}]}]
    }

    SNV_copy_counter = simulated_chromosomes_to_SNV_counts(simulated_chromosomes)

    expected_SNV_copy_counter = {
        "A": {"1": 1, "2": 1},
        "B": {"1": 1, "2": 1, "3": 1, "4": 1},
        "C": {"1": 1, "2": 2, "3": 1}
    }

    assert SNV_copy_counter == expected_SNV_copy_counter



test_simulated_chromosomes_to_SNV_counts()


def SNV_counts_to_SNV_multiplicities(SNV_copy_counter):
    """
    Convert the copy number counts of each SNV into SNV multiplicities.

    :param SNV_copy_counter: A dictionary with the count of copy numbers for each SNV.
    :return: A dictionary with the SNV multiplicities.
    """

    multiplicities = {}

    for chrom_number in SNV_copy_counter:
        if chrom_number not in multiplicities:
            multiplicities[chrom_number] = {}

        for SNV in SNV_copy_counter[chrom_number]:
            CN = SNV_copy_counter[chrom_number][SNV]

            if CN not in multiplicities[chrom_number]:
                multiplicities[chrom_number][CN] = 1
            else:
                multiplicities[chrom_number][CN] += 1

    return multiplicities


def count_SNV_multiplicities(simulated_chromosomes):
    return SNV_counts_to_SNV_multiplicities(simulated_chromosomes_to_SNV_counts(simulated_chromosomes))


# Test case for SNV_counts_to_SNV_multiplicities
def test_SNV_counts_to_SNV_multiplicities():

    test_cases = [
        {
            "input": {
                "SNV_copy_counter": {
                    "A": {"1": 1, "2": 1},
                    "B": {"1": 1, "2": 1, "3": 1, "4": 1}
                }
            },
            "expected_output": {
                "A": {1: 2},
                "B": {1: 4}
            }
        },
        {
            "input": {
                "SNV_copy_counter": {
                    "A": {"1": 2},
                    "B": {"2": 2, "3": 1, "4": 1}
                }
            },
            "expected_output": {
                "A": {2: 1},
                "B": {1: 2, 2: 1}
            }
        },
        {
            "input": {
                "SNV_copy_counter": {
                    "A": {},
                    "B": {}
                }
            },
            "expected_output": {
                "A": {},
                "B": {}
            }
        },
        {
            "input": {
                "SNV_copy_counter": {
                    "A": {"1": 1},
                    "B": {"2": 2, "3": 1, "4": 1}
                }
            },
            "expected_output": {
                "A": {1: 1},
                "B": {2: 1, 1: 2}
            }
        }
    ]

    for i, test_case in enumerate(test_cases):
        input_data = test_case["input"]
        expected_output = test_case["expected_output"]
        multiplicities = SNV_counts_to_SNV_multiplicities(input_data["SNV_copy_counter"])

        assert multiplicities == expected_output, f"Test case {i} failed: Expected {expected_output}, but got {multiplicities}."



test_SNV_counts_to_SNV_multiplicities()



##### STEP 4; now we have SNV counts, make all possible trees that could explain those SNV counts for the given epoch structure “(pre,mid, post)”
##### 
##### 
##### 
##### 
##### 
##### To save redundancy and speed up computational time-complexity of discovering and sorting through all evolutionary 
##### trees that explain the SNV data we do not insert copy number 1’s when finding all tree structures too iterate over.
##### They can be considered implicit and inserted everywhere at the very end. 
##### Non unitary copy numbers are not necessarily 
##### The reason for this is that there is a bijection between tree structures that include copy number 1’s and 
##### have paired timeline arrays where it every copy number must have existed for non zero time except for copy number 1
##### and the tree structures that do not include copy number 1 as every leaf of every tree and force all nodes to have 
##### a non zero evolutionary history where SNVs were allowed to accumulate. The latter is easier to computationally discover.
##### The tree structures can be simply constructed in a tuple of tuples format like tree = (value, left subtree, right subtree). 

### some basic tree functions:
def sort_tree(tree):
    """
    Sorts a tree in ascending order.

    Args:
        tree (tuple or None): The tree to be sorted.

    Returns:
        tuple or None: The sorted tree.

    """
    if tree is None:
        return None

    if len(tree) == 1:
        return tree

    if len(tree) == 2:
        value, subtree = tree
        return (value, sort_tree(subtree))

    if len(tree) == 3:
        value, left_subtree, right_subtree = tree
        sorted_left = sort_tree(left_subtree)
        sorted_right = sort_tree(right_subtree)

        if sorted_left[0] > sorted_right[0]:
            return (value, sorted_left, sorted_right)
        else:
            return (value, sorted_right, sorted_left)

def test_sort_tree():
    # Test case 1: Empty tree
    tree = None
    assert sort_tree(tree) == None

    # Test case 2: Tree with a single value
    tree = (5,)
    assert sort_tree(tree) == (5,)

    # Test case 3: Tree with two values
    tree = (4, (2,))
    assert sort_tree(tree) == (4, (2,))

    # Test case 4: Tree with three values
    tree = (3, (5,), (2,))
    assert sort_tree(tree) == (3, (5,), (2,))

    # Test case 5: Tree with nested values
    tree = (7, (4, (8,), (6,)), (9,))
    assert sort_tree(tree) == (7, (9,), (4, (8,), (6,)))

test_sort_tree()


def compare_trees(tree1, tree2):
    tree1 = sort_tree(tree1)
    tree2 = sort_tree(tree2)
    len_tree1 = len(tree1)
    len_tree2 = len(tree2)

    if tree1 is None and tree2 is None:
        return True
    if tree1 is None or tree2 is None:
        return False

    # If the lengths of the trees are different, they are not equal
    if len_tree1 != len_tree2:
        return False

    # If the lengths of the trees are both 1, compare their values
    if len_tree1 == 1 and len_tree2 == 1:
        return tree1[0] == tree2[0]

    # If the lengths of the trees are both 2, compare their values
    if len_tree1 == 2 and len_tree2 == 2:
        if tree1[0] != tree2[0]:
            return False
        return compare_trees(tree1[1], tree2[1])

    # If the lengths of the trees are both 3, compare their values
    if len_tree1 == 3 and len_tree2 == 3:
        if tree1[0] != tree2[0]:
            return False
        return compare_trees(tree1[1], tree2[1]) and compare_trees(tree1[2], tree2[2])

    # If the code reaches this point, the trees are not equal
    return False

def forests_are_equal(trees1, trees2):
    if len(trees1) != len(trees2):
        return False
    for tree1 in trees1:
        found_match = False
        for tree2 in trees2:
            if compare_trees(tree1, tree2):
                found_match = True
                break
        if not found_match:
            return False
    return True

def test_compare_trees_and_forests_are_equal():
    # Test cases for compare_trees
    tree1 = (5, (3,), (2, (1,), (1,)))
    tree2 = (5, (3,), (2, (1,), (1,)))
    tree3 = (5, (2, (1,), (1,)), (3,))
    tree4 = (5, (4,), (1, (1,), (0,)))
    tree5 = (5, (1, (1,), (0,)), (4,))

    assert compare_trees(tree1, tree2), "Expected tree1 and tree2 to be equal"
    assert compare_trees(tree1, tree3), "Expected tree1 and tree3 to be not equal"
    assert compare_trees(tree4, tree5), "Expected tree4 and tree5 to be not equal"

    # Test cases for forests_are_equal
    trees1 = [tree1, tree3]
    trees2 = [tree2, tree3]
    trees3 = [tree1, tree4]
    trees4 = [tree1, tree2, tree3, tree4]
    trees5 = [tree2, tree3, tree1, tree4]

    assert forests_are_equal(trees1, trees2), "Expected trees1 and trees2 to be equal"
    assert forests_are_equal(trees1, trees3), "Expected trees1 and trees3 to be not equal"
    assert forests_are_equal(trees4, trees5), "Expected trees4 and trees5 to be equal"


# Call the test function
test_compare_trees_and_forests_are_equal()


def insert_node(trees, CN):
    """
    Insert a node with the given copy number (CN) into every possible location in the input trees.

    :param trees: A list of binary trees.
    :param CN: The copy number to be inserted into the trees.
    :return: A list of new trees with the CN inserted once into each input tree.
    """

    # base case
    if trees == [] or trees == [()]:
        return [(CN,)]
    new_trees = []

    # otherwise
    for tree in trees:
        if len(tree) == 1 and CN < tree[0]:
            new_CNs = (CN, tree[0] - CN)
            new_trees.append(sort_tree((tree[0], (max(new_CNs),), (min(new_CNs),))))
        elif len(tree) == 3:
            for subtree in insert_node([tree[1]], CN):
                new_trees.append(sort_tree((tree[0], subtree, tree[2])))

            for subtree in insert_node([tree[2]], CN):
                new_trees.append(sort_tree((tree[0], tree[1], subtree)))

    return new_trees


def test_insert_node():
    def sort_expected_output(trees):
        return [sort_tree(tree) for tree in trees]

    # Test Case 1: Inserting a node into an empty tree
    trees = [()]
    CN = 5
    result = insert_node(trees, CN)
    expected_output = [(5,)]
    assert result == sort_expected_output(expected_output), f"Expected {expected_output}, but got {result}"

    # Test Case 2: Inserting a node into a tree with only one node
    trees = [(10,)]
    CN = 3
    result = insert_node(trees, CN)
    expected_output = [(10, (7,), (3,))]
    assert result == sort_expected_output(expected_output), f"Expected {expected_output}, but got {result}"

    # Test Case 3: Inserting a node into a tree with multiple nodes
    trees = [(12, (5,), (7,))]
    CN = 2
    result = insert_node(trees, CN)
    expected_output = [
                        (12, (5, (3,), (2,)), (7,)),
                        (12, (5,), (7, (5,), (2,)))
                      ]
    assert result == sort_expected_output(expected_output), f"Expected {expected_output}, but got {result}"

    trees = [(7,)]
    CN = 10
    result = insert_node(trees, CN)
    expected_output = []
    assert result == sort_expected_output(expected_output), f"Expected {expected_output}, but got {result}"

    # Test Case 6: Inserting a node with a CN that's equal to the root node
    trees = [(5,)]
    CN = 5
    result = insert_node(trees, CN)
    expected_output = []
    assert result == sort_expected_output(expected_output), f"Expected {expected_output}, but got {result}"


test_insert_node()

from general_functions import *


def complete_tree(tree):
    """
    Given a tree structure as tuple, completes the tree by adding nodes and
    making it left-heavy. Returns a set of completed tree structures.

    :param tree: A tuple representing a tree structure with 1 to 3 elements.
    :return: A set of completed tree structures.
    :raises ValueError: If the input tree tuple has a length other than 1, 2, or 3.
    """
    def _complete_tree_three_elements(tree):
        left_trees = complete_tree(tree[1])
        right_trees = complete_tree(tree[2])

        completed_trees = set()
        for left_tree in left_trees:
            for right_tree in right_trees:
                completed_trees.add(make_left_heavy((tree[0], left_tree, right_tree)))
        return completed_trees


    def _complete_tree_two_elements(tree):
        left_trees = complete_tree(tree[1])
        completed_trees = set()

        for left_tree in left_trees:
            for i in range(1, tree[0] - left_tree[0]):
                right_trees = complete_tree((tree[0] - left_tree[0] - i,))
                for right_tree in right_trees:
                    completed_trees.add(make_left_heavy((tree[0], left_tree, right_tree)))
        return completed_trees


    def _complete_tree_one_element(tree):
        if tree[0] == 0 or tree[0] == 1:
            return {tree}

        completed_trees = set()
        for i in range(1, tree[0]):
            left_trees = complete_tree((i,))
            right_trees = complete_tree((tree[0] - i,))

            for left_tree in left_trees:
                for right_tree in right_trees:
                    completed_trees.add(make_left_heavy((tree[0], left_tree, right_tree)))
        return completed_trees

    if len(tree) == 3:
        return _complete_tree_three_elements(tree)
    elif len(tree) == 2:
        return _complete_tree_two_elements(tree)
    elif len(tree) == 1:
        return _complete_tree_one_element(tree)
    else:
        raise ValueError("Invalid tree structure. Expected tuple length between 1 and 3.")


def test_complete_tree():
    tree = (4,)
    completed_trees = complete_tree(tree)
    expected_trees = {
        (4, (3, (2, (1,), (1,)), (1,)), (1,)),
        (4, (2, (1,), (1,)), (2, (1,), (1,)))
    }
    assert completed_trees == expected_trees, f"Expected: {expected_trees}, Got: {completed_trees}"

    tree = (5,)
    completed_trees = complete_tree(tree)
    expected_trees = {
        (5, (3, (2, (1,), (1,)), (1,)), (2, (1,), (1,))),
        (5, (4, (2, (1,), (1,)), (2, (1,), (1,))), (1,)),
        (5, (4, (3, (2, (1,), (1,)), (1,)), (1,)), (1,))
    }
    assert completed_trees == expected_trees, f"Expected: {expected_trees}, Got: {completed_trees}"

    tree = ((5, (3, (2,), (1,)), (2,)),)
    completed_trees = complete_tree(tree)
    expected_trees = {(5, (3, (2, (1,), (1,)), (1,)), (2, (1,), (1,)))}

    assert completed_trees == expected_trees, f"Expected: {expected_trees}, Got: {completed_trees}"

#test_complete_tree()


def complete_trees(trees):
    """
    Given a list of trees, returns a new list with each tree completed.
    
    Args:
    trees (list): A list of trees, where each tree is a tuple.
    
    Returns:
    list: A list of completed trees.
    """
    result = set()

    for tree in trees:
        result = result.union(complete_tree(tree))

    return list(result)


def test_complete_trees():
    trees = [
        (2, (1,), (1,)),
        (3, (1,), (2,)),
        (4, (2,), (2,))
    ]
    
    completed_trees_result = complete_trees(trees)
    expected_result = [
        (2, (1,), (1,)),
        (3, (2, (1,), (1,)), (1,)),
        (4, (2, (1,), (1,)), (2, (1,), (1,)))
    ]
    
    assert completed_trees_result == expected_result, f"Expected {expected_result}, but got {completed_trees_result}"

    trees = [
        (6, (4,), (2,)), 
        (6, (4, (2,), (2,)), (2,))
    ]
    completed_trees_result = complete_trees(trees)
    expected_result = [
        (6, (4, (3, (2, (1,), (1,)), (1,)), (1,)), (2, (1,), (1,))), 
        (6, (4, (2, (1,), (1,)), (2, (1,), (1,))), (2, (1,), (1,)))
    ]

    assert completed_trees_result == expected_result, f"Expected {expected_result}, but got {completed_trees_result}"
    
test_complete_trees()


def generate_trees(observed_CNs, SNV_CNs):
    """
    Generate all possible trees of every copy number from the multiplicity counts of the chromosomes and the SNVs.

    :param observed_CNs: A list of observed copy numbers.
    :param SNV_CNs: A list of copy numbers for SNVs.
    :return: A list of unique trees.
    """

    SNV_CNs.sort(reverse=True)
    observed_CNs.sort(reverse=True)

    # Initialize trees with the following tree for each chromosome
    trees = [(sum(observed_CNs), (max(observed_CNs),), (min(observed_CNs),))]

    for SNV_CN in SNV_CNs:
        if SNV_CN == 1:
            continue

        trees_with_new_node = insert_node(trees, SNV_CN)

        if not trees_with_new_node:
            assert SNV_CN in observed_CNs
            continue

        if SNV_CN in observed_CNs:
            trees += trees_with_new_node
        else:
            trees = trees_with_new_node

        while True:
            trees_with_node_inserted_again = insert_node(trees_with_new_node, SNV_CN)
            if not trees_with_node_inserted_again:
                break

            trees += trees_with_node_inserted_again
            trees_with_new_node = trees_with_node_inserted_again

    # Insert the "leaf nodes" into the tree - all of which are of CN 1
    trees = complete_trees(trees)

    trees = [sort_tree(tree) for tree in trees] 
    return trees



def test_generate_trees():
    # Test case 1
    observed_CNs_1 = [3, 2]
    SNV_CNs_1 = [3, 2]
    expected_trees_1 = [(5, (3, (2, (1,), (1,)), (1,)), (2, (1,), (1,)))]
    expected_trees_1 = [sort_tree(tree) for tree in expected_trees_1]
    result_trees_1 = generate_trees(observed_CNs_1, SNV_CNs_1)
    assert forests_are_equal(result_trees_1, expected_trees_1), f"Expected {expected_trees_1}, but got {result_trees_1}"

    # Test case 2
    observed_CNs_2 = [4, 1]
    SNV_CNs_2 = [2]
    expected_trees_2 = [(5, (4, (2, (1,), (1,)), (2, (1,), (1,))), (1,))]
    expected_trees_2 = [sort_tree(tree) for tree in expected_trees_2]
    result_trees_2 = generate_trees(observed_CNs_2, SNV_CNs_2)
    assert forests_are_equal(result_trees_2, expected_trees_2), f"Expected {expected_trees_2}, but got {result_trees_2}"

    observed_CNs_3 = [3, 0]
    SNV_CNs_3 = [3, 2, 1]
    result_trees_3 = generate_trees(observed_CNs_3, SNV_CNs_3)
    expected_trees_3 = [(3, (3, (2, (1,), (1,)), (1,)), (0,))]
    expected_trees_3 = [sort_tree(tree) for tree in expected_trees_3]
    assert forests_are_equal(result_trees_3, expected_trees_3), f"Expected {expected_trees_3}, but got {result_trees_3}"


# Call the test function
test_generate_trees()



##### STEP 5; now that all trees have been created, calculate all possible timings for each tree
##### 
##### 
##### 
##### 
##### 


def label_tree(tree, label_count, parents, label_to_copy):
    if label_to_copy == {} or label_count == 0:
        assert parents == {}, f"Unexpected value for parents: {parents}"
        assert label_to_copy == {}, f"Unexpected value for label_to_copy: {label_to_copy}"
        assert label_count == 0, f"Unexpected value for label_count: {label_count}"

    tree = list(tree)

    unique_label = label_count
    label_to_copy[unique_label] = tree[0]
    tree[0] = unique_label

    new_parent = unique_label

    if len(tree) >= 2:
        tree[1], label_count, parents, label_to_copy = label_tree(tree[1], label_count+1, parents, label_to_copy)
        parents[tree[1][0]] = new_parent

        if len(tree) == 3:
            tree[2], label_count, parents, label_to_copy = label_tree(tree[2], label_count+1, parents, label_to_copy)
            parents[tree[2][0]] = new_parent

    return (tree, label_count, parents, label_to_copy)


def test_label_tree():
    tree_1 = (5, (3,), (2, (1,), (1,)))
    expected_tree_1 = (0, (1,), (2, (3,), (4,)))
    parents_1 = {1: 0, 2: 0, 3: 2, 4: 2}
    label_to_copy_1 = {0: 5, 1: 3, 2: 2, 3: 1, 4: 1}

    result_tree_1, _, result_parents_1, result_label_to_copy_1 = label_tree(tree_1, 0, {}, {})

    assert forests_are_equal([result_tree_1], [expected_tree_1]), f"Expected {expected_tree_1}, but got {result_tree_1}"
    assert result_parents_1 == parents_1, f"Expected {parents_1}, but got {result_parents_1}"
    assert result_label_to_copy_1 == label_to_copy_1, f"Expected {label_to_copy_1}, but got {result_label_to_copy_1}"

    tree_2 = (6, (3,), (3, (2,), (1,)))
    expected_tree_2 = (0, (1,), (2, (3,), (4,)))
    parents_2 = {1: 0, 2: 0, 3: 2, 4: 2}
    label_to_copy_2 = {0: 6, 1: 3, 2: 3, 3: 2, 4: 1}

    result_tree_2, _, result_parents_2, result_label_to_copy_2 = label_tree(tree_2, 0, {}, {})

    assert forests_are_equal([result_tree_2], [expected_tree_2]), f"Expected {expected_tree_2}, but got {result_tree_2}"
    assert result_parents_2 == parents_2, f"Expected {parents_2}, but got {result_parents_2}"
    assert result_label_to_copy_2 == label_to_copy_2, f"Expected {label_to_copy_2}, but got {result_label_to_copy_2}"



test_label_tree()


def initialize_epochs_created(num_labels, root_label):
    """
    Initialize a 2D numpy array with shape (1, num_labels) containing None values. The value of the root_label-th element
    in the first row is set to 0 and all other values are set to None.

    :param num_labels: The number of labels in the tree.
    :param root_label: The label of the root node.
    :return: A 2D numpy array of None values with shape (1, num_labels).
    """
    assert isinstance(num_labels, int) and num_labels > 0, "Num labels should be a positive integer."
    assert isinstance(root_label, int) and 0 <= root_label < num_labels, "Root label should be within the range of labels."
    epochs_created = np.full((1, num_labels), None)
    epochs_created[0, root_label] = -1

    assert epochs_created.shape[0] > 0, "The resulting array must have at least one row."

    return epochs_created

def test_initialize_epochs_created():
    test_cases = [
        {"num_labels": 3, "root_label": 0, "expected": [[-1, None, None]]},
        {"num_labels": 7, "root_label": 0, "expected": [[-1, None, None, None, None, None, None]]},
        {"num_labels": 9, "root_label": 0, "expected": [[-1, None, None, None, None, None, None, None, None]]},
        {"num_labels": 11, "root_label": 0, "expected": [[-1, None, None, None, None, None, None, None, None, None, None]]}
    ]

    for test_case in test_cases:
        num_labels = test_case["num_labels"]
        root_label = test_case["root_label"]
        expected = test_case["expected"]
        
        timings = initialize_epochs_created(num_labels, root_label)
        assert timings.shape == (1, num_labels), f"Expected (1, {num_labels}), but got {timings.shape}"
        assert np.all(timings == expected), f"Expected {expected}, but got {timings}"

test_initialize_epochs_created()


def handle_other_nodes(epochs_created, label_to_copy, label, parent, total_epochs):
    assert epochs_created.ndim == 2, "Timings should be a 2-dimensional array."
    assert 0 <= label < epochs_created.shape[1], "Label should be within the range of timings."
    #assert isinstance(total_epochs, int) and total_epochs >= 0, "Epochs should be a positive integer."
    assert isinstance(total_epochs, (int, np.integer)) and total_epochs >= 0, "Epochs should be a positive integer or a non-negative numpy integer."
    assert epochs_created.shape[0] > 0, "The epochs_created array must have at least one row."

    new_epochs_created = epochs_created
    for row in range(len(epochs_created)):
        parents_time = epochs_created[row][parent]
        assert(parents_time is not None)
        logging.debug(f'row:{row}')
        logging.debug(f'total_epochs:{total_epochs}')
        logging.debug(f'parent:{parent}')
        logging.debug(f'parents_time:{parents_time}')
        logging.debug(f'total_epochs:{total_epochs}')
        logging.debug(f'copynumber:{label_to_copy[label]}')

        # first handle the two nodes under the root:
        if parent == 0:
            epochs_created_temp = np.tile(epochs_created[row], (1, 1))
            epochs_created_temp[:, label] = 0
        else:    
            if parents_time <= total_epochs and label_to_copy[label] == 1:
                if parents_time != -1:
                    epochs_created_temp = np.tile(epochs_created[row], (total_epochs - parents_time + 1, 1))
                    epochs_created_temp[:, label] = list(range(parents_time, total_epochs + 1))
                else:
                    epochs_created_temp = np.tile(epochs_created[row], (total_epochs - parents_time, 1))
                    epochs_created_temp[:, label] = list(range(parents_time + 1, total_epochs + 1))
                logging.debug(f'epochs_created_temp:{epochs_created_temp}')

            elif parents_time < total_epochs:
                epochs_created_temp = np.tile(epochs_created[row], (total_epochs - parents_time, 1))
                epochs_created_temp[:, label] = list(range(parents_time + 1, total_epochs + 1))
                logging.debug(f'epochs_created_temp:{epochs_created_temp}')

            else:
                continue

        if row == 0:
            new_epochs_created = epochs_created_temp
            logging.debug(f':new_epochs_created{new_epochs_created}')
        else:
            new_epochs_created = np.vstack([new_epochs_created,epochs_created_temp])

    return new_epochs_created


def test_handle_other_nodes():
    test_cases = [
        {
            "input": {
                "epochs_created": np.array([[-1, 0, None, None, None, None, None, None, None],
                                            [-1, 1, None, None, None, None, None, None, None],
                                            [-1, 2, None, None, None, None, None, None, None]], dtype=object),
                "label_to_copy": {0: 4, 1: 4, 2: 2, 3: 1, 4: 1, 5: 2, 6: 1, 7: 1, 8: 0},
                "label": 2,
                "parent": 1,
                "max_epoch": 2
            },
            "expected": np.array([[-1, 0, 1, None, None, None, None, None, None],
                                  [-1, 0, 2, None, None, None, None, None, None],
                                  [-1, 1, 2, None, None, None, None, None, None]], dtype=object)
        },
        {
            "input": {
                "epochs_created": np.array([[-1, 0, 1, None, None, None, None, None, None],
                                            [-1, 0, 2, None, None, None, None, None, None],
                                            [-1, 1, 2, None, None, None, None, None, None]], dtype=object),
                "label_to_copy": {0: 4, 1: 4, 2: 2, 3: 1, 4: 1, 5: 2, 6: 1, 7: 1, 8: 0},
                "label": 3,
                "parent": 2,
                "max_epoch": 2
            },
            "expected": np.array([[-1, 0, 1, 1, None, None, None, None, None],
                                  [-1, 0, 1, 2, None, None, None, None, None],
                                  [-1, 0, 2, 2, None, None, None, None, None],
                                  [-1, 1, 2, 2, None, None, None, None, None]], dtype=object)
        },
        {
            "input": {
                "epochs_created": np.array([[-1, 0, 1, 1, 1, 1, None, None, None],
                                            [-1, 0, 1, 2, 2, 1, None, None, None],
                                            [-1, 0, 2, 2, 2, 2, None, None, None],
                                            [-1, 1, 2, 2, 2, 2, None, None, None]], dtype=object),
                "label_to_copy": {0: 4, 1: 4, 2: 2, 3: 1, 4: 1, 5: 2, 6: 1, 7: 1, 8: 0},
                "label": 6,
                "parent": 5,
                "max_epoch": 2
            },
            "expected": np.array([[-1, 0, 1, 1, 1, 1, 1, None, None],
                                  [-1, 0, 1, 1, 1, 1, 2, None, None],
                                  [-1, 0, 1, 2, 2, 1, 1, None, None],
                                  [-1, 0, 1, 2, 2, 1, 2, None, None],
                                  [-1, 0, 2, 2, 2, 2, 2, None, None],
                                  [-1, 1, 2, 2, 2, 2, 2, None, None]], dtype=object)
            }
        ]

    for i, test_case in enumerate(test_cases):
        result = handle_other_nodes(test_case["input"]["epochs_created"],
                                    test_case["input"]["label_to_copy"],
                                    test_case["input"]["label"],
                                    test_case["input"]["parent"],
                                    test_case["input"]["max_epoch"])

        assert np.array_equal(result, test_case["expected"]), f"Test case {i+1} failed: expected {test_case['expected']}, but got {result}"

test_handle_other_nodes()

import numpy as np

def group_columns_by_parent(parents):
    """
    Group columns by parent.

    Args:
        parents (dict): A dictionary with keys as child indices and values as their parent indices.

    Returns:
        dict: A dictionary with keys as parent indices and values as lists of their child indices.
    """
    grouped_columns = {}
    for child, parent in parents.items():
        if parent in grouped_columns:
            grouped_columns[parent].append(child)
        else:
            grouped_columns[parent] = [child]
    return grouped_columns


def test_group_columns_by_parent():
    test_cases = [
        {
            "input": {'parents': {3: 2, 4: 2, 2: 1, 5: 1, 1: 0, 6: 0}},
            "expected": {2: [3, 4], 1: [2, 5], 0: [1, 6]}
        },
        {
            "input": {'parents': {2: 1, 3: 1, 1: 0, 5: 4, 6: 4, 4: 0}},
            "expected": {1: [2, 3], 0: [1, 4], 4: [5, 6]}
        },
        {
            "input": {'parents': {3: 2, 4: 2, 2: 1, 5: 1, 1: 0, 6: 0}},
            "expected": {2: [3, 4], 1: [2, 5], 0: [1, 6]}
        },
        {
            "input": {'parents': {3: 2, 4: 2, 2: 1, 5: 1, 1: 0, 7: 6, 8: 6, 6: 0}},
            "expected": {2: [3, 4], 1: [2, 5], 0: [1, 6], 6: [7, 8]}
        },
        {
            "input": {'parents': {4: 3, 5: 3, 3: 2, 6: 2, 2: 1, 7: 1, 1: 0, 8: 0}},
            "expected": {3: [4, 5], 2: [3, 6], 1: [2, 7], 0: [1, 8]}
        }
    ]

    for i, test_case in enumerate(test_cases):
        result = group_columns_by_parent(test_case["input"]["parents"])

        assert result == test_case["expected"], f"Test case {i+1} failed: expected {test_case['expected']}, but got {result}"


test_group_columns_by_parent()




def is_valid_row(row, grouped_columns):
    """
    Check if a row is valid based on the constraints of grouped_columns.

    Args:
        row (np.array): A NumPy array representing a row of data.
        grouped_columns (dict): A dictionary with keys as parent indices and values as lists of their child indices.

    Returns:
        bool: True if the row is valid, False otherwise.
    """
    for parent, children in grouped_columns.items():
        child_columns = [row[child] for child in children]
        if len(set(child_columns)) > 1:
            return False
    return True


# this function may be redundant now...
def filter_rows_based_on_parents(epochs_created, parents):
    """
    Filter rows in the timings array based on constraints from the parents dictionary.

    Args:
        timings (np.array): A 2D NumPy array representing the input data.
        parents (dict): A dictionary with keys as child indices and values as their parent indices.

    Returns:
        np.array: A filtered 2D NumPy array containing only the valid rows based on the constraints.
    """
    grouped_columns = group_columns_by_parent(parents)
    filtered_rows = [row for row in epochs_created if is_valid_row(row, grouped_columns)]
    filtered_timings = np.array(filtered_rows, dtype=object)
    return filtered_timings

def test_filter_rows_based_on_parents():
    epochs_created = np.array([[-1, 0, 0],
                        [-1, 0, 1],
                        [-1, 1, 0],
                        [-1, 1, 1]], dtype=object)

    parents = {1: 0, 2: 0}

    expected_filtered_epochs_created = np.array([[-1, 0, 0],
                                          [-1, 1, 1]], dtype=object)

    filtered_epochs_created = filter_rows_based_on_parents(epochs_created = epochs_created, parents = parents)
    assert np.array_equal(filtered_epochs_created, expected_filtered_epochs_created), f"Expected {expected_filtered_epochs_created}, but got {filtered_epochs_created}"

# Run the test function
test_filter_rows_based_on_parents()


def find_sibling(parents: dict, query_key: int) -> list:
    """
    Given a dictionary of 'parents' containing key-value pairs where the key is a child node and 
    the value is its parent node, this function takes a 'query_key' and returns a list of sibling keys.
    If the 'query_key' does not exist in the dictionary or there are no siblings, an empty list is returned.

    :param parents: A dictionary where keys are child nodes and values are parent nodes.
    :param query_key: The key for which we are searching sibling keys.
    :return: A list of sibling keys.
    """

    # Get the parent value for the query_key.
    query_parent = parents.get(query_key)

    if query_parent is None:
        return []

    # Find sibling keys by iterating through the parents dictionary.
    sibling_keys = [
        key for key, parent in parents.items() if parent == query_parent and key != query_key
    ]
    assert(len(sibling_keys) == 1)
    return sibling_keys[0]


def test_find_sibling():
    parents = {1: 0, 2: 0, 4: 3, 5: 3}
    key = 1
    siblings = find_sibling(parents, key)
    expected_output = 2
    assert siblings == expected_output, f"Expected {expected_output}, but got {siblings}"

    key = 4
    siblings = find_sibling(parents, key)
    expected_output = 5
    assert siblings == expected_output, f"Expected {expected_output}, but got {siblings}"

test_find_sibling()


import concurrent.futures
import copy
import logging

def get_timings_per_tree(tree, total_epochs): 
    logging.debug("get_timings_per_tree")

    labelled_tree, label_count, parents, label_to_copy = label_tree(
        tree=copy.deepcopy(tree),
        label_count=0,
        parents={},
        label_to_copy={}
    )

    if total_epochs == -1:
        epochs_created = initialize_epochs_created(num_labels=label_count + 1, root_label=label)
        epochs_created = -1
        return (tree, labelled_tree, label_count, epochs_created, parents)

    for label in range(label_count + 1):
        if label == 0:
            epochs_created = initialize_epochs_created(num_labels=label_count + 1, root_label=label)
            assert epochs_created.shape[0] > 0, "The epochs_created array must have at least one row."
        else:
            logging.debug("labelled_tree")
            logging.debug(labelled_tree)
            logging.debug("copy numbers")
            logging.debug([label_to_copy[x] for x in range(label_count+1)])
            logging.debug("label")
            logging.debug(label)
            logging.debug("epochs_created")
            logging.debug(epochs_created)
            assert epochs_created.shape[0] > 0, "The epochs_created array must have at least one row."
            sibling = find_sibling(parents,label)
            if sibling < label:
                logging.debug("insertion by copying sibling")
                epochs_created[:,label] = epochs_created[:,sibling]
            else:
                logging.debug("insertion by finding all possible epochs")
                epochs_created = handle_other_nodes(epochs_created=epochs_created,
                                                    label_to_copy=label_to_copy,
                                                    label=label,
                                                    parent=parents[label],
                                                    total_epochs=total_epochs
                                                    )
            logging.debug("epochs_created")
            logging.debug(epochs_created)

            if None in epochs_created[:,label]:
                return (None, None, None, None, None)

    logging.debug("epochs_created")
    logging.debug(epochs_created)
    logging.debug("epochs_created after filter_rows")
    epochs_created = filter_rows_based_on_parents(epochs_created=epochs_created, parents=parents)
    logging.debug(epochs_created)

    return (tree, labelled_tree, label_count, epochs_created, parents)

def time_limited_get_timings_per_tree(tree, total_epochs, max_time=1.0):
    # Create a process pool executor
    with concurrent.futures.ProcessPoolExecutor() as executor:
        # Run get_timings_per_tree with the specified timeout
        future = executor.submit(get_timings_per_tree, tree, total_epochs)
        try:
            result = future.result(timeout=max_time)
        except concurrent.futures.TimeoutError:
            result = (None, None, None, None, None)
    return result



def max_tree_depth(tree):
    """
    Calculate the maximum depth of a given tree.

    The tree is represented as a tuple of nested tuples, where each inner tuple
    represents a subtree. The depth of the tree is the length of the longest
    path from the root to a leaf node.

    Args:
        tree (tuple): The tree to calculate the depth of.

    Returns:
        int: The maximum depth of the tree.
    """
    if not tree:
        return 0

    if not isinstance(tree, tuple):
        raise ValueError("Invalid input. Tree must be a tuple of nested tuples.")

    max_depth = 0
    for subtree in tree[1:]:
        depth = max_tree_depth(subtree)
        max_depth = max(max_depth, depth)

    return max_depth + 1


def test_max_tree_depth():
    tree1 = (6, (3,), (3, (2,), (1,)))
    tree2 = (1, (2, (3,)), (4, (5, (6,))))
    tree3 = (1,)
    tree4 = ()
    tree5 = (1, (2, (3, (4, (5, (6,))))))
    tree6 = (4, (2, (1,), (1,)), (2, (1,), (1,)))

    assert max_tree_depth(tree1) == 3, f"Failed for tree1: {tree1}"
    assert max_tree_depth(tree2) == 4, f"Failed for tree2: {tree2}"
    assert max_tree_depth(tree3) == 1, f"Failed for tree3: {tree3}"
    assert max_tree_depth(tree4) == 0, f"Failed for tree4: {tree4}"
    assert max_tree_depth(tree5) == 6, f"Failed for tree5: {tree5}"
    assert max_tree_depth(tree6) == 3, f"Failed for tree5: {tree5}"


test_max_tree_depth()

def tree_in_bounds(tree,total_epochs,tree_flexibility):
    # it is plus three because the root node is fictional, and the next two nodes have occured prior to the simulation and the leaf nodes don't need any SNVs
    # ahh it should actually be 2 because the root node is fictional but either the next two nodes or the root node something needs to happen
    depth = max_tree_depth(tree)
    return depth >= total_epochs + 2 - tree_flexibility and depth <= total_epochs + 2


def get_all_trees(observed_SNV_multiplicities, observed_CNs, total_epochs, tree_flexibility):
    trees = {}
    for chrom in observed_SNV_multiplicities:
        logging.debug(chrom)
        logging.debug(observed_CNs)
        logging.debug(observed_SNV_multiplicities)

        all_trees = generate_trees(
            observed_CNs=observed_CNs[chrom],
            SNV_CNs=list(observed_SNV_multiplicities[chrom].keys())
        )
        trees_to_keep = []
        for tree in all_trees:
            depth = max_tree_depth(tree)
            if tree_in_bounds(tree,total_epochs,tree_flexibility):
                trees_to_keep.append(tree)
        trees[chrom] = copy.deepcopy(trees_to_keep)

        # check there is a tree and explain why:
        if len(trees_to_keep) == 0:
            results = all_trees
            logging.info("%s", results)
            for result in results:
                logging.debug("%s", max_tree_depth(result))

            logging.info("trees are of length 0")
            tests = [tree_in_bounds(tree,total_epochs,tree_flexibility) for tree in all_trees]
            logging.info("tests %s", tests)
            logging.info("any tests %s", any(tests))
            logging.info("chrom %s", chrom)
            logging.info("trees %s", all_trees)
            logging.info("max_tree_depth %s", [max_tree_depth(tree) for tree in all_trees])
            logging.info("bounds %s", (total_epochs +2  - tree_flexibility,total_epochs + 2))
            logging.info("total_epochs %s", total_epochs)

            return None

    return trees

def get_all_timings(trees, observed_SNV_multiplicities, total_epochs, pre, mid, post):
    trees_and_timings = {}
    for chrom in observed_SNV_multiplicities:
        logging.debug(f"(pre,mid,post) + {(pre,mid,post)}")
        trees_and_timings[chrom] = get_chrom_trees_and_timings(all_trees = trees[chrom], total_epochs=total_epochs)
        logging.debug(trees_and_timings[chrom])

    return trees_and_timings



def test_get_all_trees_and_timings():
    # Fill in your input data here
    input_data = {
        'observed_SNV_multiplicities': {
            0: {3: 62, 1: 27, 2: 22},
            1: {5: 24, 3: 28, 1: 208, 2: 85},
            2: {5: 45, 2: 65, 1: 192, 3: 21},
            3: {3: 42, 2: 83, 1: 74},
            4: {2: 36, 1: 34},
            5: {4: 27, 2: 55, 1: 88, 3: 37},
            6: {2: 39, 1: 35},
            7: {1: 57, 3: 41, 2: 23},
            8: {4: 39, 2: 62, 1: 37},
            9: {1: 82},
            10: {3: 29, 1: 23, 2: 18},
            11: {2: 65, 1: 80},
            12: {3: 20, 1: 42, 2: 16},
            13: {3: 22, 1: 38, 2: 30},
            14: {4: 13, 3: 14, 1: 38, 2: 17},
            15: {1: 28},
            16: {3: 11, 2: 9, 1: 11},
            17: {4: 6, 2: 19, 1: 62},
            18: {5: 7, 3: 7, 2: 13, 1: 16},
            19: {2: 10, 1: 13},
            20: {2: 16, 1: 10, 3: 9},
            21: {3: 4, 2: 4, 1: 20},
            22: {4: 21, 2: 42, 1: 72}
        },
        'observed_CNs': {
            0: [3, 0],
            1: [5, 2],
            2: [5, 1],
            3: [3, 2],
            4: [2, 0],
            5: [4, 3],
            6: [2, 0],
            7: [3, 1],
            8: [4, 2],
            9: [1, 1],
            10: [3, 0],
            11: [2, 2],
            12: [3, 1],
            13: [3, 2],
            14: [4, 0],
            15: [1, 0],
            16: [3, 0],
            17: [4, 1],
            18: [5, 0],
            19: [2, 0],
            20: [3, 2],
            21: [3, 0],
            22: [4, 0]
        },
        'pre': 1,
        'mid': 1,
        'post': -1,
        'max_epoch': 2
    }
 

    # Call the function with the input data
    result = get_all_trees_and_timings(**input_data)

    # Fill in your expected output here
    expected_output = {0: [((3, (3, (2, (1,), (1,)), (1,)), (0,)), [0, [1, [2, [3], [4]], [5]], [6]], 6, array([[-1, 0, 1, 1, 1, 1, 0],
       [-1, 0, 1, 2, 2, 1, 0],
       [-1, 0, 2, 2, 2, 2, 0],
       [-1, 1, 2, 2, 2, 2, 1]], dtype=object), {3: 2, 4: 2, 2: 1, 5: 1, 1: 0, 6: 0})], 1: [((7, (5, (3, (2, (1,), (1,)), (1,)), (2, (1,), (1,))), (2, (1,), (1,))), [0, [1, [2, [3, [4], [5]], [6]], [7, [8], [9]]], [10, [11], [12]]], 12, array([[-1, 0, 1, 2, 2, 2, 2, 1, 1, 1, 0, 0, 0],
       [-1, 0, 1, 2, 2, 2, 2, 1, 1, 1, 0, 1, 1],
       [-1, 0, 1, 2, 2, 2, 2, 1, 1, 1, 0, 2, 2],
       [-1, 0, 1, 2, 2, 2, 2, 1, 2, 2, 0, 0, 0],
       [-1, 0, 1, 2, 2, 2, 2, 1, 2, 2, 0, 1, 1],
       [-1, 0, 1, 2, 2, 2, 2, 1, 2, 2, 0, 2, 2]], dtype=object), {4: 3, 5: 3, 3: 2, 6: 2, 2: 1, 8: 7, 9: 7, 7: 1, 1: 0, 11: 10, 12: 10, 10: 0})], 2: [((6, (5, (3, (2, (1,), (1,)), (1,)), (2, (1,), (1,))), (1,)), [0, [1, [2, [3, [4], [5]], [6]], [7, [8], [9]]], [10]], 10, array([[-1, 0, 1, 2, 2, 2, 2, 1, 1, 1, 0],
       [-1, 0, 1, 2, 2, 2, 2, 1, 2, 2, 0]], dtype=object), {4: 3, 5: 3, 3: 2, 6: 2, 2: 1, 8: 7, 9: 7, 7: 1, 1: 0, 10: 0})], 3: [((5, (3, (2, (1,), (1,)), (1,)), (2, (1,), (1,))), [0, [1, [2, [3], [4]], [5]], [6, [7], [8]]], 8, array([[-1, 0, 1, 1, 1, 1, 0, 0, 0],
       [-1, 0, 1, 1, 1, 1, 0, 1, 1],
       [-1, 0, 1, 1, 1, 1, 0, 2, 2],
       [-1, 0, 1, 2, 2, 1, 0, 0, 0],
       [-1, 0, 1, 2, 2, 1, 0, 1, 1],
       [-1, 0, 1, 2, 2, 1, 0, 2, 2],
       [-1, 0, 2, 2, 2, 2, 0, 0, 0],
       [-1, 0, 2, 2, 2, 2, 0, 1, 1],
       [-1, 0, 2, 2, 2, 2, 0, 2, 2],
       [-1, 1, 2, 2, 2, 2, 1, 1, 1],
       [-1, 1, 2, 2, 2, 2, 1, 2, 2]], dtype=object), {3: 2, 4: 2, 2: 1, 5: 1, 1: 0, 7: 6, 8: 6, 6: 0})], 4: [((2, (2, (1,), (1,)), (0,)), [0, [1, [2], [3]], [4]], 4, array([[-1, 0, 0, 0, 0],
       [-1, 0, 1, 1, 0],
       [-1, 0, 2, 2, 0],
       [-1, 1, 1, 1, 1],
       [-1, 1, 2, 2, 1],
       [-1, 2, 2, 2, 2]], dtype=object), {2: 1, 3: 1, 1: 0, 4: 0})], 5: [((7, (4, (2, (1,), (1,)), (2, (1,), (1,))), (3, (2, (1,), (1,)), (1,))), [0, [1, [2, [3], [4]], [5, [6], [7]]], [8, [9, [10], [11]], [12]]], 12, array([[-1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1],
       [-1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 2, 2, 1],
       [-1, 0, 1, 1, 1, 1, 1, 1, 0, 2, 2, 2, 2],
       [-1, 0, 1, 1, 1, 1, 2, 2, 0, 1, 1, 1, 1],
       [-1, 0, 1, 1, 1, 1, 2, 2, 0, 1, 2, 2, 1],
       [-1, 0, 1, 1, 1, 1, 2, 2, 0, 2, 2, 2, 2],
       [-1, 0, 1, 2, 2, 1, 1, 1, 0, 1, 1, 1, 1],
       [-1, 0, 1, 2, 2, 1, 1, 1, 0, 1, 2, 2, 1],
       [-1, 0, 1, 2, 2, 1, 1, 1, 0, 2, 2, 2, 2],
       [-1, 0, 1, 2, 2, 1, 2, 2, 0, 1, 1, 1, 1],
       [-1, 0, 1, 2, 2, 1, 2, 2, 0, 1, 2, 2, 1],
       [-1, 0, 1, 2, 2, 1, 2, 2, 0, 2, 2, 2, 2],
       [-1, 0, 2, 2, 2, 2, 2, 2, 0, 1, 1, 1, 1],
       [-1, 0, 2, 2, 2, 2, 2, 2, 0, 1, 2, 2, 1],
       [-1, 0, 2, 2, 2, 2, 2, 2, 0, 2, 2, 2, 2],
       [-1, 1, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 2]], dtype=object), {3: 2, 4: 2, 2: 1, 6: 5, 7: 5, 5: 1, 1: 0, 10: 9, 11: 9, 9: 8, 12: 8, 8: 0}), ((7, (4, (3, (2, (1,), (1,)), (1,)), (1,)), (3, (2, (1,), (1,)), (1,))), [0, [1, [2, [3, [4], [5]], [6]], [7]], [8, [9, [10], [11]], [12]]], 12, array([[-1, 0, 1, 2, 2, 2, 2, 1, 0, 1, 1, 1, 1],
       [-1, 0, 1, 2, 2, 2, 2, 1, 0, 1, 2, 2, 1],
       [-1, 0, 1, 2, 2, 2, 2, 1, 0, 2, 2, 2, 2]], dtype=object), {4: 3, 5: 3, 3: 2, 6: 2, 2: 1, 7: 1, 1: 0, 10: 9, 11: 9, 9: 8, 12: 8, 8: 0})], 6: [((2, (2, (1,), (1,)), (0,)), [0, [1, [2], [3]], [4]], 4, array([[-1, 0, 0, 0, 0],
       [-1, 0, 1, 1, 0],
       [-1, 0, 2, 2, 0],
       [-1, 1, 1, 1, 1],
       [-1, 1, 2, 2, 1],
       [-1, 2, 2, 2, 2]], dtype=object), {2: 1, 3: 1, 1: 0, 4: 0})], 7: [((4, (3, (2, (1,), (1,)), (1,)), (1,)), [0, [1, [2, [3], [4]], [5]], [6]], 6, array([[-1, 0, 1, 1, 1, 1, 0],
       [-1, 0, 1, 2, 2, 1, 0],
       [-1, 0, 2, 2, 2, 2, 0],
       [-1, 1, 2, 2, 2, 2, 1]], dtype=object), {3: 2, 4: 2, 2: 1, 5: 1, 1: 0, 6: 0})], 8: [((6, (4, (3, (2, (1,), (1,)), (1,)), (1,)), (2, (1,), (1,))), [0, [1, [2, [3, [4], [5]], [6]], [7]], [8, [9], [10]]], 10, array([[-1, 0, 1, 2, 2, 2, 2, 1, 0, 0, 0],
       [-1, 0, 1, 2, 2, 2, 2, 1, 0, 1, 1],
       [-1, 0, 1, 2, 2, 2, 2, 1, 0, 2, 2]], dtype=object), {4: 3, 5: 3, 3: 2, 6: 2, 2: 1, 7: 1, 1: 0, 9: 8, 10: 8, 8: 0}), ((6, (4, (2, (1,), (1,)), (2, (1,), (1,))), (2, (1,), (1,))), [0, [1, [2, [3], [4]], [5, [6], [7]]], [8, [9], [10]]], 10, array([[-1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0],
       [-1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1],
       [-1, 0, 1, 1, 1, 1, 1, 1, 0, 2, 2],
       [-1, 0, 1, 1, 1, 1, 2, 2, 0, 0, 0],
       [-1, 0, 1, 1, 1, 1, 2, 2, 0, 1, 1],
       [-1, 0, 1, 1, 1, 1, 2, 2, 0, 2, 2],
       [-1, 0, 1, 2, 2, 1, 1, 1, 0, 0, 0],
       [-1, 0, 1, 2, 2, 1, 1, 1, 0, 1, 1],
       [-1, 0, 1, 2, 2, 1, 1, 1, 0, 2, 2],
       [-1, 0, 1, 2, 2, 1, 2, 2, 0, 0, 0],
       [-1, 0, 1, 2, 2, 1, 2, 2, 0, 1, 1],
       [-1, 0, 1, 2, 2, 1, 2, 2, 0, 2, 2],
       [-1, 0, 2, 2, 2, 2, 2, 2, 0, 0, 0],
       [-1, 0, 2, 2, 2, 2, 2, 2, 0, 1, 1],
       [-1, 0, 2, 2, 2, 2, 2, 2, 0, 2, 2],
       [-1, 1, 2, 2, 2, 2, 2, 2, 1, 1, 1],
       [-1, 1, 2, 2, 2, 2, 2, 2, 1, 2, 2]], dtype=object), {3: 2, 4: 2, 2: 1, 6: 5, 7: 5, 5: 1, 1: 0, 9: 8, 10: 8, 8: 0})], 9: [((2, (1,), (1,)), [0, [1], [2]], 2, array([[-1, 0, 0],
       [-1, 1, 1],
       [-1, 2, 2]], dtype=object), {1: 0, 2: 0})], 10: [((3, (3, (2, (1,), (1,)), (1,)), (0,)), [0, [1, [2, [3], [4]], [5]], [6]], 6, array([[-1, 0, 1, 1, 1, 1, 0],
       [-1, 0, 1, 2, 2, 1, 0],
       [-1, 0, 2, 2, 2, 2, 0],
       [-1, 1, 2, 2, 2, 2, 1]], dtype=object), {3: 2, 4: 2, 2: 1, 5: 1, 1: 0, 6: 0})], 11: [((4, (2, (1,), (1,)), (2, (1,), (1,))), [0, [1, [2], [3]], [4, [5], [6]]], 6, array([[-1, 0, 0, 0, 0, 0, 0],
       [-1, 0, 0, 0, 0, 1, 1],
       [-1, 0, 0, 0, 0, 2, 2],
       [-1, 0, 1, 1, 0, 0, 0],
       [-1, 0, 1, 1, 0, 1, 1],
       [-1, 0, 1, 1, 0, 2, 2],
       [-1, 0, 2, 2, 0, 0, 0],
       [-1, 0, 2, 2, 0, 1, 1],
       [-1, 0, 2, 2, 0, 2, 2],
       [-1, 1, 1, 1, 1, 1, 1],
       [-1, 1, 1, 1, 1, 2, 2],
       [-1, 1, 2, 2, 1, 1, 1],
       [-1, 1, 2, 2, 1, 2, 2],
       [-1, 2, 2, 2, 2, 2, 2]], dtype=object), {2: 1, 3: 1, 1: 0, 5: 4, 6: 4, 4: 0})], 12: [((4, (3, (2, (1,), (1,)), (1,)), (1,)), [0, [1, [2, [3], [4]], [5]], [6]], 6, array([[-1, 0, 1, 1, 1, 1, 0],
       [-1, 0, 1, 2, 2, 1, 0],
       [-1, 0, 2, 2, 2, 2, 0],
       [-1, 1, 2, 2, 2, 2, 1]], dtype=object), {3: 2, 4: 2, 2: 1, 5: 1, 1: 0, 6: 0})], 13: [((5, (3, (2, (1,), (1,)), (1,)), (2, (1,), (1,))), [0, [1, [2, [3], [4]], [5]], [6, [7], [8]]], 8, array([[-1, 0, 1, 1, 1, 1, 0, 0, 0],
       [-1, 0, 1, 1, 1, 1, 0, 1, 1],
       [-1, 0, 1, 1, 1, 1, 0, 2, 2],
       [-1, 0, 1, 2, 2, 1, 0, 0, 0],
       [-1, 0, 1, 2, 2, 1, 0, 1, 1],
       [-1, 0, 1, 2, 2, 1, 0, 2, 2],
       [-1, 0, 2, 2, 2, 2, 0, 0, 0],
       [-1, 0, 2, 2, 2, 2, 0, 1, 1],
       [-1, 0, 2, 2, 2, 2, 0, 2, 2],
       [-1, 1, 2, 2, 2, 2, 1, 1, 1],
       [-1, 1, 2, 2, 2, 2, 1, 2, 2]], dtype=object), {3: 2, 4: 2, 2: 1, 5: 1, 1: 0, 7: 6, 8: 6, 6: 0})], 14: [((4, (4, (3, (2, (1,), (1,)), (1,)), (1,)), (0,)), [0, [1, [2, [3, [4], [5]], [6]], [7]], [8]], 8, array([[-1, 0, 1, 2, 2, 2, 2, 1, 0]], dtype=object), {4: 3, 5: 3, 3: 2, 6: 2, 2: 1, 7: 1, 1: 0, 8: 0})], 15: [((1, (1,), (0,)), [0, [1], [2]], 2, array([[-1, 0, 0],
       [-1, 1, 1],
       [-1, 2, 2]], dtype=object), {1: 0, 2: 0})], 16: [((3, (3, (2, (1,), (1,)), (1,)), (0,)), [0, [1, [2, [3], [4]], [5]], [6]], 6, array([[-1, 0, 1, 1, 1, 1, 0],
       [-1, 0, 1, 2, 2, 1, 0],
       [-1, 0, 2, 2, 2, 2, 0],
       [-1, 1, 2, 2, 2, 2, 1]], dtype=object), {3: 2, 4: 2, 2: 1, 5: 1, 1: 0, 6: 0})], 17: [((5, (4, (2, (1,), (1,)), (2, (1,), (1,))), (1,)), [0, [1, [2, [3], [4]], [5, [6], [7]]], [8]], 8, array([[-1, 0, 1, 1, 1, 1, 1, 1, 0],
       [-1, 0, 1, 1, 1, 1, 2, 2, 0],
       [-1, 0, 1, 2, 2, 1, 1, 1, 0],
       [-1, 0, 1, 2, 2, 1, 2, 2, 0],
       [-1, 0, 2, 2, 2, 2, 2, 2, 0],
       [-1, 1, 2, 2, 2, 2, 2, 2, 1]], dtype=object), {3: 2, 4: 2, 2: 1, 6: 5, 7: 5, 5: 1, 1: 0, 8: 0})], 18: [((5, (5, (3, (2, (1,), (1,)), (1,)), (2, (1,), (1,))), (0,)), [0, [1, [2, [3, [4], [5]], [6]], [7, [8], [9]]], [10]], 10, array([[-1, 0, 1, 2, 2, 2, 2, 1, 1, 1, 0],
       [-1, 0, 1, 2, 2, 2, 2, 1, 2, 2, 0]], dtype=object), {4: 3, 5: 3, 3: 2, 6: 2, 2: 1, 8: 7, 9: 7, 7: 1, 1: 0, 10: 0})], 19: [((2, (2, (1,), (1,)), (0,)), [0, [1, [2], [3]], [4]], 4, array([[-1, 0, 0, 0, 0],
       [-1, 0, 1, 1, 0],
       [-1, 0, 2, 2, 0],
       [-1, 1, 1, 1, 1],
       [-1, 1, 2, 2, 1],
       [-1, 2, 2, 2, 2]], dtype=object), {2: 1, 3: 1, 1: 0, 4: 0})], 20: [((5, (3, (2, (1,), (1,)), (1,)), (2, (1,), (1,))), [0, [1, [2, [3], [4]], [5]], [6, [7], [8]]], 8, array([[-1, 0, 1, 1, 1, 1, 0, 0, 0],
       [-1, 0, 1, 1, 1, 1, 0, 1, 1],
       [-1, 0, 1, 1, 1, 1, 0, 2, 2],
       [-1, 0, 1, 2, 2, 1, 0, 0, 0],
       [-1, 0, 1, 2, 2, 1, 0, 1, 1],
       [-1, 0, 1, 2, 2, 1, 0, 2, 2],
       [-1, 0, 2, 2, 2, 2, 0, 0, 0],
       [-1, 0, 2, 2, 2, 2, 0, 1, 1],
       [-1, 0, 2, 2, 2, 2, 0, 2, 2],
       [-1, 1, 2, 2, 2, 2, 1, 1, 1],
       [-1, 1, 2, 2, 2, 2, 1, 2, 2]], dtype=object), {3: 2, 4: 2, 2: 1, 5: 1, 1: 0, 7: 6, 8: 6, 6: 0})], 21: [((3, (3, (2, (1,), (1,)), (1,)), (0,)), [0, [1, [2, [3], [4]], [5]], [6]], 6, array([[-1, 0, 1, 1, 1, 1, 0],
       [-1, 0, 1, 2, 2, 1, 0],
       [-1, 0, 2, 2, 2, 2, 0],
       [-1, 1, 2, 2, 2, 2, 1]], dtype=object), {3: 2, 4: 2, 2: 1, 5: 1, 1: 0, 6: 0})], 22: [((4, (4, (2, (1,), (1,)), (2, (1,), (1,))), (0,)), [0, [1, [2, [3], [4]], [5, [6], [7]]], [8]], 8, array([[-1, 0, 1, 1, 1, 1, 1, 1, 0],
       [-1, 0, 1, 1, 1, 1, 2, 2, 0],
       [-1, 0, 1, 2, 2, 1, 1, 1, 0],
       [-1, 0, 1, 2, 2, 1, 2, 2, 0],
       [-1, 0, 2, 2, 2, 2, 2, 2, 0],
       [-1, 1, 2, 2, 2, 2, 2, 2, 1]], dtype=object), {3: 2, 4: 2, 2: 1, 6: 5, 7: 5, 5: 1, 1: 0, 8: 0})]} 

    # Assert that the function output matches the expected output
    assert result == expected_output

def test_get_timings_per_tree():
    test_cases = [
        {
            'input': {'tree': (6, (4, (2, (1,), (1,)), (2, (1,), (1,))), (2, (1,), (1,))), 'max_epoch': 2},
            'output': ((6, (4, (2, (1,), (1,)), (2, (1,), (1,))), (2, (1,), (1,))), [0, [1, [2, [3], [4]], [5, [6], [7]]], [8, [9], [10]]], 10, np.array([[-1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0],
       [-1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1],
       [-1, 0, 1, 1, 1, 1, 1, 1, 0, 2, 2],
       [-1, 0, 1, 1, 1, 1, 2, 2, 0, 0, 0],
       [-1, 0, 1, 1, 1, 1, 2, 2, 0, 1, 1],
       [-1, 0, 1, 1, 1, 1, 2, 2, 0, 2, 2],
       [-1, 0, 1, 2, 2, 1, 1, 1, 0, 0, 0],
       [-1, 0, 1, 2, 2, 1, 1, 1, 0, 1, 1],
       [-1, 0, 1, 2, 2, 1, 1, 1, 0, 2, 2],
       [-1, 0, 1, 2, 2, 1, 2, 2, 0, 0, 0],
       [-1, 0, 1, 2, 2, 1, 2, 2, 0, 1, 1],
       [-1, 0, 1, 2, 2, 1, 2, 2, 0, 2, 2],
       [-1, 0, 2, 2, 2, 2, 2, 2, 0, 0, 0],
       [-1, 0, 2, 2, 2, 2, 2, 2, 0, 1, 1],
       [-1, 0, 2, 2, 2, 2, 2, 2, 0, 2, 2],
       [-1, 1, 2, 2, 2, 2, 2, 2, 1, 1, 1],
       [-1, 1, 2, 2, 2, 2, 2, 2, 1, 2, 2]], dtype=object), {3: 2, 4: 2, 2: 1, 6: 5, 7: 5, 5: 1, 1: 0, 9: 8, 10: 8, 8: 0})
        },
        {
            'input': {'tree': (3, (3, (2, (1,), (1,)), (1,)), (0,)), 'max_epoch': 2},
            'output': ((3, (3, (2, (1,), (1,)), (1,)), (0,)), [0, [1, [2, [3], [4]], [5]], [6]], 6, np.array([[-1, 0, 1, 1, 1, 1, 0],
       [-1, 0, 1, 2, 2, 1, 0],
       [-1, 0, 2, 2, 2, 2, 0],
       [-1, 1, 2, 2, 2, 2, 1]], dtype=object), {3: 2, 4: 2, 2: 1, 5: 1, 1: 0, 6: 0})
        },
        {
            'input': {'tree': (5, (3, (2, (1,), (1,)), (1,)), (2, (1,), (1,))), 'max_epoch': 2},
            'output': ((5, (3, (2, (1,), (1,)), (1,)), (2, (1,), (1,))), [0, [1, [2, [3], [4]], [5]], [6, [7], [8]]], 8, np.array([[-1, 0, 1, 1, 1, 1, 0, 0, 0],
       [-1, 0, 1, 1, 1, 1, 0, 1, 1],
       [-1, 0, 1, 1, 1, 1, 0, 2, 2],
       [-1, 0, 1, 2, 2, 1, 0, 0, 0],
       [-1, 0, 1, 2, 2, 1, 0, 1, 1],
       [-1, 0, 1, 2, 2, 1, 0, 2, 2],
       [-1, 0, 2, 2, 2, 2, 0, 0, 0],
       [-1, 0, 2, 2, 2, 2, 0, 1, 1],
       [-1, 0, 2, 2, 2, 2, 0, 2, 2],
       [-1, 1, 2, 2, 2, 2, 1, 1, 1],
       [-1, 1, 2, 2, 2, 2, 1, 2, 2]], dtype=object), {3: 2, 4: 2, 2: 1, 5: 1, 1: 0, 7: 6, 8: 6, 6: 0})
        },
        {
            'input': {'tree': (3, (3, (2, (1,), (1,)), (1,)), (0,)), 'max_epoch': 2},
            'output': ((3, (3, (2, (1,), (1,)), (1,)), (0,)), [0, [1, [2, [3], [4]], [5]], [6]], 6, np.array([[-1, 0, 1, 1, 1, 1, 0],
       [-1, 0, 1, 2, 2, 1, 0],
       [-1, 0, 2, 2, 2, 2, 0],
       [-1, 1, 2, 2, 2, 2, 1]], dtype=object), {3: 2, 4: 2, 2: 1, 5: 1, 1: 0, 6: 0})
        }
    ]

    for test_case in test_cases:
        input_data = test_case['input']
        expected_output = test_case['output']
        result = get_timings_per_tree(**input_data)
        assert result == expected_output, f'For {input_data}, expected {expected_output} but got {result}'

# Call the test function
#test_get_timings_per_tree()




def get_chrom_trees_and_timings(all_trees,total_epochs):
    logging.debug("all trees")
    logging.debug(all_trees)

    #chrom_trees_and_timings = [get_timings_per_tree(tree=x, total_epochs=total_epochs) for x in all_trees]
    #chrom_trees_and_timings = [time_limited_get_timings_per_tree(tree=x, total_epochs=total_epochs) for x in all_trees]
    #get_timings_per_tree, tree, total_epochs

    chrom_trees_and_timings = [get_timings_per_tree(tree=x, total_epochs=total_epochs) for x in all_trees]

    logging.debug("chrom_trees_and_timings")
    logging.debug(chrom_trees_and_timings)

    chrom_trees_and_timings = [x for x in chrom_trees_and_timings if x[3] is not None and not None in x[3]]

    logging.debug("chrom_trees_and_timings")
    logging.debug(chrom_trees_and_timings)

    #assert(len(chrom_trees_and_timings) != 0)

    return chrom_trees_and_timings

def test_get_chrom_trees_and_timings():
    test_cases = [
        {
            'input': {'tree': (6, (4, (2, (1,), (1,)), (2, (1,), (1,))), (2, (1,), (1,))), 'max_epoch': 2},
            'output': ((6, (4, (2, (1,), (1,)), (2, (1,), (1,))), (2, (1,), (1,))), [0, [1, [2, [3], [4]], [5, [6], [7]]], [8, [9], [10]]], 10, np.array([[-1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0],
       [-1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1],
       [-1, 0, 1, 1, 1, 1, 1, 1, 0, 2, 2],
       [-1, 0, 1, 1, 1, 1, 2, 2, 0, 0, 0],
       [-1, 0, 1, 1, 1, 1, 2, 2, 0, 1, 1],
       [-1, 0, 1, 1, 1, 1, 2, 2, 0, 2, 2],
       [-1, 0, 1, 2, 2, 1, 1, 1, 0, 0, 0],
       [-1, 0, 1, 2, 2, 1, 1, 1, 0, 1, 1],
       [-1, 0, 1, 2, 2, 1, 1, 1, 0, 2, 2],
       [-1, 0, 1, 2, 2, 1, 2, 2, 0, 0, 0],
       [-1, 0, 1, 2, 2, 1, 2, 2, 0, 1, 1],
       [-1, 0, 1, 2, 2, 1, 2, 2, 0, 2, 2],
       [-1, 0, 2, 2, 2, 2, 2, 2, 0, 0, 0],
       [-1, 0, 2, 2, 2, 2, 2, 2, 0, 1, 1],
       [-1, 0, 2, 2, 2, 2, 2, 2, 0, 2, 2],
       [-1, 1, 2, 2, 2, 2, 2, 2, 1, 1, 1],
       [-1, 1, 2, 2, 2, 2, 2, 2, 1, 2, 2]], dtype=object), {3: 2, 4: 2, 2: 1, 6: 5, 7: 5, 5: 1, 1: 0, 9: 8, 10: 8, 8: 0})
        },
        {
            'input': {'tree': (3, (3, (2, (1,), (1,)), (1,)), (0,)), 'max_epoch': 2},
            'output': ((3, (3, (2, (1,), (1,)), (1,)), (0,)), [0, [1, [2, [3], [4]], [5]], [6]], 6, np.array([[-1, 0, 1, 1, 1, 1, 0],
       [-1, 0, 1, 2, 2, 1, 0],
       [-1, 0, 2, 2, 2, 2, 0],
       [-1, 1, 2, 2, 2, 2, 1]], dtype=object), {3: 2, 4: 2, 2: 1, 5: 1, 1: 0, 6: 0})
        },
        {
            'input': {'tree': (5, (3, (2, (1,), (1,)), (1,)), (2, (1,), (1,))), 'max_epoch': 2},
            'output': ((5, (3, (2, (1,), (1,)), (1,)), (2, (1,), (1,))), [0, [1, [2, [3], [4]], [5]], [6, [7], [8]]], 8, np.array([[-1, 0, 1, 1, 1, 1, 0, 0, 0],
       [-1, 0, 1, 1, 1, 1, 0, 1, 1],
       [-1, 0, 1, 1, 1, 1, 0, 2, 2],
       [-1, 0, 1, 2, 2, 1, 0, 0, 0],
       [-1, 0, 1, 2, 2, 1, 0, 1, 1],
       [-1, 0, 1, 2, 2, 1, 0, 2, 2],
       [-1, 0, 2, 2, 2, 2, 0, 0, 0],
       [-1, 0, 2, 2, 2, 2, 0, 1, 1],
       [-1, 0, 2, 2, 2, 2, 0, 2, 2],
       [-1, 1, 2, 2, 2, 2, 1, 1, 1],
       [-1, 1, 2, 2, 2, 2, 1, 2, 2]], dtype=object), {3: 2, 4: 2, 2: 1, 5: 1, 1: 0, 7: 6, 8: 6, 6: 0})
        },
        {
            'input': {'tree': (3, (3, (2, (1,), (1,)), (1,)), (0,)), 'max_epoch': 2},
            'output': ((3, (3, (2, (1,), (1,)), (1,)), (0,)), [0, [1, [2, [3], [4]], [5]], [6]], 6, np.array([[-1, 0, 1, 1, 1, 1, 0],
       [-1, 0, 1, 2, 2, 1, 0],
       [-1, 0, 2, 2, 2, 2, 0],
       [-1, 1, 2, 2, 2, 2, 1]], dtype=object), {3: 2, 4: 2, 2: 1, 5: 1, 1: 0, 6: 0})
        }
    ]

    for test_case in test_cases:
        input_data = test_case['input']
        expected_output = test_case['output']
        result = get_chrom_trees_and_timings(**input_data)
        assert result == expected_output, f'For {input_data}, expected {expected_output} but got {result}'

# Call the test function
#test_get_chrom_trees_and_timings()


# Test cases for get_all_trees_and_timings function
def test_get_all_trees_and_timings():
    observed_SNV_multiplicities = {"chr1": {1: 1, 2: 1}}
    observed_CNs = {"chr1": [1, 2]}
    pre, mid, post = 1, 1, 1
    max_epoch = calculate_epochs(pre,mid,post) 
    trees_and_timings = get_all_trees_and_timings(observed_SNV_multiplicities, observed_CNs, pre, mid, post,max_epoch)

    assert "chr1" in trees_and_timings, "Chromosome key is missing in the result"
    assert len(trees_and_timings["chr1"]) > 0, "No trees and timings found for the chromosome"


def find_indices(list_to_check, item_to_find):
    indices = [i for i, x in enumerate(list_to_check) if x == item_to_find]
    return indices

def test_find_indices():
    test_list = [1, 2, 3, 2, 4, 2, 5]
    item = 2
    indices = find_indices(test_list, item)
    assert indices == [1, 3, 5], "find_indices returned incorrect indices"


test_find_indices()


def find_non_parent_children(parents):
    """
    Find children that are not parents in the given parent-child relationship dictionary.

    Args:
        parents (dict): A dictionary with keys as child indices and values as their parent indices.

    Returns:
        set: A set of child indices that are not parents.
    """
    all_children = set(parents.keys())
    all_parents = set(parents.values())
    non_parent_children = all_children.difference(all_parents)
    return non_parent_children

def test_find_non_parent_children():
    parents = {1: 0, 2: 0, 4: 3, 5: 3}
    expected_non_parent_children = {1, 2, 4, 5}
    non_parent_children = find_non_parent_children(parents)
    assert non_parent_children == expected_non_parent_children, f"Expected {expected_non_parent_children}, but got {non_parent_children}"

# Run the test function
test_find_non_parent_children()


def calculate_child_parent_diff(epochs_created: np.ndarray, parents: Dict, max_epoch: int):
    """
    Calculate the difference between child and parent epochs.

    Args:
        epochs_created (np.ndarray): A 2D numpy array containing the epochs when elements were created.
        parents (Dict): A dictionary with keys as child indices and values as their parent indices.
        epochs (int): The total number of epochs.

    Returns:
        np.ndarray: A 2D numpy array containing the difference between child and parent epochs.
    """
    assert isinstance(epochs_created, np.ndarray) and epochs_created.ndim == 2, f"epochs_created must be a 2D numpy array, but was {type(epochs_created)} with dimensions {epochs_created.ndim}. epochs_created: {str(epochs_created)}"

    assert isinstance(parents, dict), f"parents must be a dictionary, but was {type(parents)}"
    
    branch_lengths = np.copy(epochs_created)
    for child in parents: 
        branch_lengths[:, parents[child]] = epochs_created[:, child] - epochs_created[:, parents[child]]

    non_parent_children = find_non_parent_children(parents)
    for column in non_parent_children:
        branch_lengths[:, column] = max_epoch - epochs_created[:, column]
        
    return branch_lengths

def test_calculate_child_parent_diff():
    # First test case
    epochs_created = np.array([[-1, 0, 1, 1, 1, 1, 0],
                               [-1, 0, 1, 2, 2, 1, 0],
                               [-1, 0, 2, 2, 2, 2, 0],
                               [-1, 1, 2, 2, 2, 2, 1]], dtype=object)
    parents = {3: 2, 4: 2, 2: 1, 5: 1, 1: 0, 6: 0}
    max_epoch = 2
    result = calculate_child_parent_diff(epochs_created, parents, max_epoch)
    expected = np.array([[1, 1, 0, 1, 1, 1, 2],
                         [1, 1, 1, 0, 0, 1, 2],
                         [1, 2, 0, 0, 0, 0, 2],
                         [2, 1, 0, 0, 0, 0, 1]], dtype=object)
    assert np.array_equal(result, expected), f"calculate_child_parent_diff returned incorrect values.\nExpected: {expected}\nGot: {result}"

    # Second test case
    epochs_created = np.array([[-1, 0, 1, 1, 1, 1, 1, 1, 0],
                               [-1, 0, 1, 1, 1, 1, 2, 2, 0],
                               [-1, 0, 1, 2, 2, 1, 1, 1, 0],
                               [-1, 0, 1, 2, 2, 1, 2, 2, 0],
                               [-1, 0, 2, 2, 2, 2, 2, 2, 0],
                               [-1, 1, 2, 2, 2, 2, 2, 2, 1]], dtype=object)
    parents = {3: 2, 4: 2, 2: 1, 6: 5, 7: 5, 5: 1, 1: 0, 8: 0}
    max_epoch = 2
    result = calculate_child_parent_diff(epochs_created, parents, max_epoch)
    expected = np.array([[1, 1, 0, 1, 1, 0, 1, 1, 2],
                         [1, 1, 0, 1, 1, 1, 0, 0, 2],
                         [1, 1, 1, 0, 0, 0, 1, 1, 2],
                         [1, 1, 1, 0, 0, 1, 0, 0, 2],
                         [1, 2, 0, 0, 0, 0, 0, 0, 2],
                         [2, 1, 0, 0, 0, 0, 0, 0, 1]], dtype=object)
    assert np.array_equal(result, expected), f"calculate_child_parent_diff returned incorrect values.\nExpected: {expected}\nGot: {result}"


test_calculate_child_parent_diff()

def extract_copy_numbers(tree):
    CNs = [x for x in re.split("\(|\)|,|'", str(tree)) if x.isdigit()]
    CNs = [int(x) for x in CNs]
    return CNs


def test_extract_copy_numbers():
    tree = "(3,(1,1))"
    CNs = extract_copy_numbers(tree)
    assert CNs == [3, 1, 1], "extract_copy_numbers returned incorrect CNs"


test_extract_copy_numbers()


def stack_same_CN_branch_lengths(unique_CNs, CNs, branch_lengths):
    for i, CN in enumerate(unique_CNs):
        indices = find_indices(CNs, CN)
        new_stacked_branch_lengths = branch_lengths[:, indices].sum(axis=1)

        if i == 0:
            stacked_branch_lengths = new_stacked_branch_lengths
        else:
            stacked_branch_lengths = np.vstack((stacked_branch_lengths, new_stacked_branch_lengths))

    return np.transpose(stacked_branch_lengths)

def stack_same_CN_branch_lengths(CNs, branch_lengths):
    """
    Stacks branch lengths of the same copy number state.

    Arguments:
    unique_CNs -- a list of unique copy number states.
    CNs -- a list of copy number states, matching the second dimension of branch_lengths.
    branch_lengths -- a 2D array where the second dimension matches CNs.

    Returns:
    A 2D numpy array where branch lengths of the same copy number state have been stacked.
    """
    # Remove the first column
    CNs = CNs[1:]
    branch_lengths = branch_lengths[:, 1:]
    unique_CNs = sorted(list(set(CNs)), reverse = True)

    for i, CN in enumerate(unique_CNs):
        indices = find_indices(CNs, CN)
        new_stacked_branch_lengths = branch_lengths[:, indices].sum(axis=1)

        if i == 0:
            stacked_branch_lengths = new_stacked_branch_lengths
        else:
            stacked_branch_lengths = np.vstack((stacked_branch_lengths, new_stacked_branch_lengths))

    return np.transpose(stacked_branch_lengths), unique_CNs


def stack_same_CN_branch_lengths(CNs, branch_lengths):
    """
    Stacks branch lengths of the same copy number state.

    Arguments:
    CNs -- a list of copy number states, matching the second dimension of branch_lengths.
    branch_lengths -- a 2D array where the second dimension matches CNs.

    Returns:
    A 2D numpy array where branch lengths of the same copy number state have been stacked.
    """
    # Remove the first column
    CNs = CNs[1:]
    branch_lengths = branch_lengths[:, 1:]
    unique_CNs = sorted(list(set(CNs)), reverse=True)

    stacked_branch_lengths = []  # Initialize an empty list

    for CN in unique_CNs:
        indices = find_indices(CNs, CN)
        new_stacked_branch_lengths = branch_lengths[:, indices].sum(axis=1)
        stacked_branch_lengths.append(new_stacked_branch_lengths)

    stacked_branch_lengths = np.vstack(stacked_branch_lengths).T

    return stacked_branch_lengths, unique_CNs



def test_stack_same_CN_branch_lengths():
    # First test case
    unique_CNs = [5, 3, 2, 1]
    CNs = [5, 3, 2, 1, 1, 1, 2, 1, 1]
    branch_lengths = np.array([[1, 1, 0, 1, 1, 1, 0, 2, 2],
                               [1, 1, 0, 1, 1, 1, 1, 1, 1],
                               [1, 1, 0, 1, 1, 1, 2, 0, 0],
                               [1, 1, 1, 0, 0, 1, 0, 2, 2],
                               [1, 1, 1, 0, 0, 1, 1, 1, 1],
                               [1, 1, 1, 0, 0, 1, 2, 0, 0],
                               [1, 2, 0, 0, 0, 0, 0, 2, 2],
                               [1, 2, 0, 0, 0, 0, 1, 1, 1],
                               [1, 2, 0, 0, 0, 0, 2, 0, 0],
                               [2, 1, 0, 0, 0, 0, 0, 1, 1],
                               [2, 1, 0, 0, 0, 0, 1, 0, 0]], dtype=object)
    result = stack_same_CN_branch_lengths(CNs, branch_lengths)
    expected = np.array([[1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2],
                         [1, 1, 1, 1, 1, 1, 2, 2, 2, 1, 1],
                         [0, 1, 2, 1, 2, 3, 0, 1, 2, 0, 1],
                         [7, 5, 3, 5, 3, 1, 4, 2, 0, 2, 0]], dtype=object).T
    assert np.array_equal(result, expected), f"stack_same_CN_branch_lengths returned incorrect values.\nExpected: {expected}\nGot: {result}"

    # Second test case
    unique_CNs = [3, 2, 1, 0]
    CNs = [3, 3, 2, 1, 1, 1, 0]
    branch_lengths = np.array([[1, 1, 0, 1, 1, 1, 2],
                               [1, 1, 1, 0, 0, 1, 2],
                               [1, 2, 0, 0, 0, 0, 2],
                               [2, 1, 0, 0, 0, 0, 1]], dtype=object)
    result = stack_same_CN_branch_lengths(CNs, branch_lengths)
   
    expected = np.array([[2, 2, 3, 3],
                         [0, 1, 0, 0],
                         [3, 1, 0, 0],
                         [2, 2, 2, 1]], dtype=object).T
    assert np.array_equal(result, expected), f"stack_same_CN_branch_lengths returned incorrect values.\nExpected: {expected}\nGot: {result}"

    # Third test case
    unique_CNs = [4, 2, 1, 0]
    CNs = [4, 4, 2, 1, 1, 2, 1, 1, 0]
    branch_lengths = np.array([[1, 1, 0, 1, 1, 0, 1, 1, 2],
                               [1, 1, 0, 1, 1, 1, 0, 0, 2],
                               [1, 1, 1, 0, 0, 0, 1, 1, 2],
                               [1, 1, 1, 0, 0, 1, 0, 0, 2],
                               [1, 2, 0, 0, 0, 0, 0, 0, 2],
                               [2, 1, 0, 0, 0, 0, 0, 0, 1]], dtype=object)
    result = stack_same_CN_branch_lengths(CNs, branch_lengths)
    expected = np.array([[2, 2, 2, 2, 3, 3],
                         [0, 1, 1, 2, 0, 0],
                         [4, 2, 2, 0, 0, 0],
                         [2, 2, 2, 2, 2, 1]], dtype=object).T
    assert np.array_equal(result, expected), f"stack_same_CN_branch_lengths returned incorrect values.\nExpected: {expected}\nGot: {result}"


#test_stack_same_CN_branch_lengths()



def get_branch_lengths(trees_and_timings, max_epoch):
    tree, labelled_tree, count, epochs_created, parents = trees_and_timings
    branch_lengths = calculate_child_parent_diff(
        epochs_created = epochs_created, 
        parents = parents,
        max_epoch = max_epoch
    )
    CNs = extract_copy_numbers(tree)
    stacked_branch_lengths,unique_CNs = stack_same_CN_branch_lengths(CNs, branch_lengths)
   
    logging.debug("trees_and_timings")
    logging.debug(str(trees_and_timings))
    logging.debug("max_epoch")
    logging.debug(max_epoch)
    logging.debug("branch_lengths")
    logging.debug(branch_lengths)
    logging.debug("stacked_branch_lengths")
    logging.debug(stacked_branch_lengths)

    return CNs, unique_CNs, branch_lengths, stacked_branch_lengths


# Test cases for get_branch_lengths function
def test_get_branch_lengths():
    # First test case
    trees_and_timings = ((5, (4, (2, (1,), (1,)), (2, (1,), (1,))), (1,)), [0, [1, [2, [3], [4]], [5, [6], [7]]], [8]], 8, np.array([[-1, 0, 1, 1, 1, 1, 1, 1, 0],
       [-1, 0, 1, 1, 1, 1, 2, 2, 0],
       [-1, 0, 1, 2, 2, 1, 1, 1, 0],
       [-1, 0, 1, 2, 2, 1, 2, 2, 0],
       [-1, 0, 2, 2, 2, 2, 2, 2, 0],
       [-1, 1, 2, 2, 2, 2, 2, 2, 1]], dtype=object), {3: 2, 4: 2, 2: 1, 6: 5, 7: 5, 5: 1, 1: 0, 8: 0})
    max_epoch = 2
    CNs, unique_CNs, branch_lengths, stacked_branch_lengths = get_branch_lengths(trees_and_timings, max_epoch)

    assert CNs == [5, 4, 2, 1, 1, 2, 1, 1, 1]
    assert unique_CNs == [4, 2, 1]
    assert np.array_equal(branch_lengths, np.array([[1, 1, 0, 1, 1, 0, 1, 1, 2],
                                                     [1, 1, 0, 1, 1, 1, 0, 0, 2],
                                                     [1, 1, 1, 0, 0, 0, 1, 1, 2],
                                                     [1, 1, 1, 0, 0, 1, 0, 0, 2],
                                                     [1, 2, 0, 0, 0, 0, 0, 0, 2],
                                                     [2, 1, 0, 0, 0, 0, 0, 0, 1]], dtype=object))
    assert np.array_equal(stacked_branch_lengths, np.array([[1, 1, 0, 6],
                                                             [1, 1, 1, 4],
                                                             [1, 1, 1, 4],
                                                             [1, 1, 2, 2],
                                                             [1, 2, 0, 2],
                                                             [2, 1, 0, 1]], dtype=object))


    # Second test case
    trees_and_timings = ((5, (5, (3, (2, (1,), (1,)), (1,)), (2, (1,), (1,))), (0,)), [0, [1, [2, [3, [4], [5]], [6]], [7, [8], [9]]], [10]], 10, np.array([[-1, 0, 1, 2, 2, 2, 2, 1, 1, 1, 0],
       [-1, 0, 1, 2, 2, 2, 2, 1, 2, 2, 0]], dtype=object), {4: 3, 5: 3, 3: 2, 6: 2, 2: 1, 8: 7, 9: 7, 7: 1, 1: 0, 10: 0})
    max_epoch = 2
    CNs, unique_CNs, branch_lengths, stacked_branch_lengths = get_branch_lengths(trees_and_timings, max_epoch)

    assert CNs == [5, 5, 3, 2, 1, 1, 1, 2, 1, 1, 0]
    assert unique_CNs == [5, 3, 2, 1, 0]
    assert np.array_equal(branch_lengths, np.array([[1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 2],
                                                     [1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 2]], dtype=object))
    assert np.array_equal(stacked_branch_lengths, np.array([[2, 1, 0, 2, 2],
                                                             [2, 1, 1, 0, 2]], dtype=object))

    # Third test case
    trees_and_timings = ((2, (2, (1,), (1,)), (0,)), [0, [1, [2], [3]], [4]], 4, np.array([[-1, 0, 0, 0, 0],
       [-1, 0, 1, 1, 0],
       [-1, 0, 2, 2, 0],
       [-1, 1, 1, 1, 1],
       [-1, 1, 2, 2, 1],
       [-1, 2, 2, 2, 2]], dtype=object), {2: 1, 3: 1, 1: 0, 4: 0})
    max_epoch = 2
    CNs, unique_CNs, branch_lengths, stacked_branch_lengths = get_branch_lengths(trees_and_timings, max_epoch)

    assert CNs == [2, 2, 1, 1, 0]
    assert unique_CNs == [2, 1, 0]
    assert np.array_equal(branch_lengths, np.array([[1, 0, 2, 2, 2],
                                                     [1, 1, 1, 1, 2],
                                                     [1, 2, 0, 0, 2],
                                                     [2, 0, 1, 1, 1],
                                                     [2, 1, 0, 0, 1],
                                                     [3, 0, 0, 0, 0]], dtype=object))
    assert np.array_equal(stacked_branch_lengths, np.array([[1, 4, 2],
                                                             [2, 2, 2],
                                                             [3, 0, 2],
                                                             [2, 2, 1],
                                                             [3, 0, 1],
                                                             [3, 0, 0]], dtype=object))

#test_get_branch_lengths()

##### STEP 7; from the branch lengths calculate the BP likelihoods
##### 
##### 
##### 
##### 
##### 
def get_path_code(code_list):
    output = ""
    count = 0

    for code in code_list:
        if code == "A":
            count += 1
        elif code == "GD":
            output += str(count)
            count = 0
            output += "G"

    output += str(count)
    return output


def test_get_path_code():
    # Test case 1
    code_list = ["A", "A", "GD", "A", "A", "A", "GD", "A"]
    result = get_path_code(code_list)
    expected = "2G3G1"
    assert result == expected, f"get_path_code returned incorrect result: {result}"

    # Test case 2
    code_list = []
    result = get_path_code(code_list)
    expected = "0"
    assert result == expected, f"get_path_code returned incorrect result: {result}"

    # Test case 3
    code_list = ['A', 'GD']
    result = get_path_code(code_list)
    expected = '1G0'
    assert result == expected, f"get_path_code returned incorrect result: {result}"

test_get_path_code()


# Configure logging settings (you only need to do this once in your script or module)
# this would be a good idea to use throughout the script

def timing_struct_to_all_structures(trees_and_timings, pre, mid, post, max_epoch):
    all_structures = {}
    
    for chrom in trees_and_timings:
        all_structures[chrom] = timing_structs_to_all_structs_per_chrom(trees_and_timings[chrom], pre, mid, post, max_epoch)

    

    return all_structures


def timing_structs_to_all_structs_per_chrom(trees_and_timings, pre, mid, post, max_epoch):
    logging.debug("trees_and_timings")
    logging.debug(trees_and_timings)
    all_structures = [] 
    for index, these_tts in enumerate(trees_and_timings):  # these tts are a 2d array
        if None in these_tts[3]:
            continue
            #BP_likelihoods = -1
        else:
            # trace back to here, asdfasdf
            CNs, unique_CNs, branch_lengths, stacked_branch_lengths = get_branch_lengths(trees_and_timings=these_tts, max_epoch=max_epoch)

            logging.debug("CNs, unique_CNs, branch_lengths, stacked_branch_lengths")
            logging.debug(f"{CNs}, {unique_CNs}, {branch_lengths}, {stacked_branch_lengths}")
            logging.debug("starts and ends")

            path = create_path(pre, mid, post)

            logging.debug(path)

            starts = these_tts[3] #+1
            ends = these_tts[3] + branch_lengths #+1

            logging.debug("starts")
            logging.debug(starts)
            logging.debug("ends")
            logging.debug(ends)

            paths = calculate_BP_paths(branch_lengths, starts, ends, path)

        tree, labelled_tree, count, epochs_created, parents = these_tts

        all_structures += [{
            "pre": pre,
            "mid": mid,
            "post": post,
            "path": path,
            "tree": tree,
            "parents": parents,
            "labelled_tree": labelled_tree,
            "count": count,
            "epochs_created": epochs_created,
            "CNs": CNs,
            "branch_lengths": branch_lengths,
            "unique_CNs": unique_CNs,
            "stacked_branch_lengths": stacked_branch_lengths,
            "starts":starts,
            "ends":ends,
            "paths":paths
        }]

    return all_structures



def timing_struct_to_BP_likelihood_per_chrom(data, structures, p_up, p_down):
    #logging.getLogger().setLevel(logging.DEBUG)

    all_BP_likelihoods = []

    for structure in structures: 
        assert(None not in structure["epochs_created"])
        logging.debug("structures")
        logging.debug(structures)

        likelihoods = calculate_likelihoods_from_paths(paths=structure["paths"], CNs=structure["CNs"], data=data, p_up=p_up, p_down=p_down)

        logging.debug("paths")
        logging.debug(structure["paths"])

        logging.debug("likelihoods")
        logging.debug(likelihoods)

        ll = np.log(likelihoods)
        
        logging.debug("log likelihoods")
        logging.debug(ll)
        structure["BP_individual_log_likelihoods"] = ll

        BP_likelihood = np.sum(ll[:, 1:], axis=1)

        logging.debug("sum them to get BP likelihoods")
        logging.debug(BP_likelihood)

        structure["BP_loglikelihood"] = BP_likelihood

    #logging.getLogger().setLevel(logging.INFO)

def test_timing_struct_to_BP_likelihood_per_chrom():
    data = {
        "3G": [None, [0.2, 0.6, 0.2]],
        "1G1": [None, [0.9, 0.1, 0]],
        "2G0": [None, [1, 0, 0]],
    }
    trees_and_timings = [
        (
            None,
            None,
            None,
            np.array([[2, 3], [3, 2]]),
            None,
        ),
        (
            None,
            None,
            None,
            np.array([[0, 2], [2, 1]]),
            None,
        ),
    ]
    pre, mid, post = 2, 1, 1

    expected_all_BP_likelihoods = [
        np.array([-0.916290731874155, -0.916290731874155]),
        np.array([-0.40546510810816444, -0.40546510810816444]),
    ]

    all_BP_likelihoods = timing_struct_to_BP_likelihood_per_chrom(data, trees_and_timings, pre, mid, post)

    for i in range(len(all_BP_likelihoods)):
        assert np.allclose(all_BP_likelihoods[i], expected_all_BP_likelihoods[i], rtol=1e-9, atol=1e-9)


#test_timing_struct_to_BP_likelihood_per_chrom()


def create_path(pre, mid, post):
    path = []
    if pre > 0:
        path += ["A"] * pre
    if mid > -1:
        path += ["GD"]
    if mid > 0:
        path += ["A"] * mid
    if post > -1:
        path += ["GD"]
    if post > 0:
        path += ["A"] * post

    return path


def test_create_path():
    assert create_path(3, 2, 1) == ["A", "A", "A", "GD", "A", "A", "GD", "A"]
    assert create_path(0, 0, 0) == ["GD", "GD"]
    assert create_path(2, 0, 2) == ["A", "A", "GD", "GD", "A", "A"]
    assert create_path(1, 1, -1) == ["A", "GD", "A"]

test_create_path()


def calculate_BP_paths(branch_lengths, starts, ends, path):
    paths = np.zeros(ends.shape, dtype=object, order="C")

    for row in range(branch_lengths.shape[0]):
        for col in range(branch_lengths.shape[1]):

            these_paths = path[starts[row][col] : ends[row][col]]  # MODIFIED THIS, BUT NOT ENTIRELY SURE, CHECK, ERROR
            path_code = get_path_code(these_paths)

            paths[row][col] = path_code

    # this needs to be updated for paths with length longer than 1. 
    # suppose a path was like GG1. The branching process can take this and calculate a probability, however when the SNV likelihood doe sthe same thing it gives an incompatible probability. 
    # the SNV likelihood needs to know how long SNV can accumulate on a branch. for SNVs to accurately take iunto consideration of how lon g they can accumulate, the BP lieklihood needs to recognise this. 

    # therefore GG1 really needs to look like p(GG1) = p(GD, GD) * p(D)^3 * 4 * p(L)
    # this is not straightforward to program. How do we do it?
    # furthermroe, how do we recognise it?
    #

    # first we recoginise this phenomenum as occuring to paths greater than length 1. all other paths naturally work.
    # P(GG1 1 to 2 true) = P(GG1 1 to 1) * UA/(1-U-D)
    # this works because there si always one chromosome where nothing happens to it. in the 1 to 1 caseA

    # does this also work fro paths of length 1?
    # yes if it is a non genome doubling one. 
    # actually it always works unless the last epoch is a gd one. then it is just the usual. /calculate_BP_paths

    return paths

import numpy as np

def test_calculate_BP_paths():
    # Test case 1
    branch_lengths = np.array([[1, 2, 2], [2, 1, 1], [3, 0, 0]], dtype=object)
    starts = np.array([[-1, 0, 0], [-1, 1, 1], [-1, 2, 2]], dtype=object)
    ends = np.array([[0, 2, 2], [1, 2, 2], [2, 2, 2]], dtype=object)
    path = ['A', 'GD', 'A']
    result = calculate_BP_paths(branch_lengths, starts, ends, path)
    expected = np.array([['0', '1G0', '1G0'], ['0', '0G0', '0G0'], ['0', '0', '0']], dtype=object)
    assert np.array_equal(result, expected), f"calculate_BP_paths returned incorrect result: {result}"

    # Test case 2
    # Test case 2
    branch_lengths = np.array([[1, 1, 0, 1, 1, 1, 2],
                               [1, 1, 1, 0, 0, 1, 2],
                               [1, 2, 0, 0, 0, 0, 2],
                               [2, 1, 0, 0, 0, 0, 1]], dtype=object)
    starts = np.array([[-1, 0, 1, 1, 1, 1, 0],
                       [-1, 0, 1, 2, 2, 1, 0],
                       [-1, 0, 2, 2, 2, 2, 0],
                       [-1, 1, 2, 2, 2, 2, 1]], dtype=object)
    ends = np.array([[0, 1, 1, 2, 2, 2, 2],
                     [0, 1, 2, 2, 2, 2, 2],
                     [0, 2, 2, 2, 2, 2, 2],
                     [1, 2, 2, 2, 2, 2, 2]], dtype=object)
    path = ['A', 'GD', 'A']
    result = calculate_BP_paths(branch_lengths, starts, ends, path)
    expected = np.array([['0', '1', '0', '0G0', '0G0', '0G0', '1G0'],
                         ['0', '1', '0G0', '0', '0', '0G0', '1G0'],
                         ['0', '1G0', '0', '0', '0', '0', '1G0'],
                         ['0', '0G0', '0', '0', '0', '0', '0G0']], dtype=object)
    assert np.array_equal(result, expected), f"calculate_BP_paths returned incorrect result: {result}"

    # Test case 3
    branch_lengths = np.array([[1, 1, 0, 1, 1, 0, 1, 1, 2],
                               [1, 1, 0, 1, 1, 1, 0, 0, 2],
                               [1, 1, 1, 0, 0, 0, 1, 1, 2],
                               [1, 1, 1, 0, 0, 1, 0, 0, 2],
                               [1, 2, 0, 0, 0, 0, 0, 0, 2],
                               [2, 1, 0, 0, 0, 0, 0, 0, 1]], dtype=object)
    starts = np.array([[-1, 0, 1, 1, 1, 1, 1, 1, 0],
                       [-1, 0, 1, 1, 1, 1, 2, 2, 0],
                       [-1, 0, 1, 2, 2, 1, 1, 1, 0],
                       [-1, 0, 1, 2, 2, 1, 2, 2, 0],
                       [-1, 0, 2, 2, 2, 2, 2, 2, 0],
                       [-1, 1, 2, 2, 2, 2, 2, 2, 1]], dtype=object)
    ends = np.array([[0, 1, 1, 2, 2, 1, 2, 2, 2],
                     [0, 1, 1, 2, 2, 2, 2, 2, 2],
                     [0, 1, 2, 2, 2, 1, 2, 2, 2],
                     [0, 1, 2, 2, 2, 2, 2, 2, 2],
                     [0, 2, 2, 2, 2, 2, 2, 2, 2],
                     [1, 2, 2, 2, 2, 2, 2, 2, 2]], dtype=object)
    path = ['A', 'GD', 'A']
    result = calculate_BP_paths(branch_lengths, starts, ends, path)
    expected = np.array([['0', '1', '0', '0G0', '0G0', '0', '0G0', '0G0', '1G0'],
                         ['0', '1', '0', '0G0', '0G0', '0G0', '0', '0', '1G0'],
                         ['0', '1', '0G0', '0', '0', '0', '0G0', '0G0', '1G0'],
                         ['0', '1', '0G0', '0', '0', '0G0', '0', '0', '1G0'],
                         ['0', '1G0', '0', '0', '0', '0', '0', '0', '1G0'],
                         ['0', '0G0', '0', '0', '0', '0', '0', '0', '0G0']], dtype=object)
    assert np.array_equal(result, expected), f"calculate_BP_paths returned incorrect result: {result}"

test_calculate_BP_paths()


def check_last_character(code):
    """
    Checks if the last character of the given code is 'G' or '0' versus a non-zero integer.

    Parameters:
        code (str): Code represented as a string.

    Returns:
        bool: True if the last character of the code is 'G' or '0', False if it's a non-zero integer.
    """

    while code and code[-1] == '0':
        code = code[:-1]
        if len(code) == 0: 
            return True

    if code and code[-1] in {'G', '0'}:
        return True
    else:
        return False


def test_check_last_character():
    """
    Test the check_last_character function with various test cases.
    """
    test_cases = {
        '0': True,
        '1': False,
        '0G0': True,
        '0G1': False,
        '1G0': True,
        '0G': True,
        '1G': True,
        'G': True,
        'G1': False,
        '10': False, 
        '01': False,
        '00': True
    }

    for code, expected in test_cases.items():
        result = check_last_character(code)
        assert result == expected, f"For code '{code}', expected {expected} but got {result}"

test_check_last_character()

def shorten_by_one(path):
    """
    Shortens the given path by one. If the last character is '0', it is removed along with the preceding 'G'. 
    If the path is '0', a ValueError is raised.
    """
    if path == '0':
        raise ValueError("Cannot shorten '0'")
    
    if path.endswith('G0'):
        return path[:-2]
    elif path[-1].isdigit() and int(path[-1]) > 0:
        return path[:-1] + str(int(path[-1]) - 1)
    else:
        raise ValueError(f"Invalid path: {path}")


def test_shorten_by_one():
    """
    Tests the `shorten_by_one` function with various path values.
    """
    test_cases = [
        ('1G2', '1G1'),
        ('1G1', '1G0'),
        ('1G0', '1'),
        ('3G4G0', '3G4'),
        ('3G4', '3G3'),
        ('3G3', '3G2'),
        ('3G2', '3G1'),
        ('3G1', '3G0'),
        ('3G0', '3'),
        ('3', '2'),
        ('2', '1'),
        ('1', '0'),
        # Testing error cases
        ('0', ValueError),
    ]

    for i, (path, expected) in enumerate(test_cases):
        try:
            result = shorten_by_one(path)
            assert result == expected, f"Test case {i+1} failed: ({path}) => {result}, expected {expected}"
        except ValueError as e:
            assert isinstance(expected, type) and issubclass(expected, Exception), f"Test case {i+1} failed: Exception was not expected for path {path}"
        except Exception as e:
            print(f"Test case {i+1} failed with an unexpected exception: {type(e).__name__} - {e}")



# Run the test function
test_shorten_by_one()


def calculate_likelihoods_from_paths(paths, CNs, data, p_up, p_down):
    likelihoods = np.zeros(paths.shape, dtype=float, order="C")

    assert isinstance(p_up, int) and 0 <= p_up <= 100, f"p_up should be an integer from 0 to 100: {p_up}"
    assert isinstance(p_down, int) and 0 <= p_down <= 100, f"p_down should be an integer from 0 to 100: {p_down}"
    p_up = p_up/100.0
    p_down = p_down/100.0

    assert 0 <= p_up <= 1 and 0 <= p_down <= 1, f"p_up: {p_up}, p_down: {p_down}, p_up and p_down should be between 0 and 1"


    for row in range(paths.shape[0]):
        for col in range(paths.shape[1]):
            path = paths[row][col]
            if check_last_character(path) or CNs[col] < 2:
                likelihood = data[path][1][min(2,CNs[col])] # the prob of going from CN 1 to ...
                # the reason for the min(2,*) is that if it is a copy number zero node then find the prob it is zero, if it is 1 then likewise, if 2 or more then it is just the prob of a bifuctation\
            else:
                # for the sake of the SNV likelihood we need the doubling up to be at the last second:
                if len(path) > 1:
                    likelihood = data[shorten_by_one(path)][1][1] * p_up
                    n = data[shorten_by_one(path)].shape[1]  # assuming 'data' is a numpy array

                    for j in range(2, n):
                        likelihood += data[shorten_by_one(path)][1][j] * p_up * (p_down ** (j-1)) * j

                else:
                    assert(len(path) == 1)
                    likelihood = data[path][1][2]



            if math.isnan(likelihood):
                logging.getLogger().setLevel(logging.DEBUG)
                logging.debug("Last element in path: %s", path[-1])
                logging.debug("Value of CNs[col]: %s", CNs[col])
                logging.debug("Value of p_up: %s", p_up)
                logging.debug("Value of p_down: %s", p_down)
                logging.debug("Calculated likelihood: %s", likelihood)
                logging.debug("Exiting the program.")
                logging.debug("Likelihood contains 'nan'. Exiting the program.")
                sys.exit()

            likelihoods[row][col] = likelihood


    if np.isnan(likelihoods).any():
        print("The likelihoods array contains NaN values. Exiting the program.")
        sys.exit()

    return likelihoods


def test_calculate_paths_and_likelihoods():
    data = {
        "3G": [None, [0.2, 0.6, 0.2]],
        "1G1": [None, [0.9, 0.1, 0]],
        "2G0": [None, [1, 0, 0]],
    }
    branch_lengths = np.array([[2, 1], [1, 1]])
    CNs = ["1", "0"]
    starts = np.array([[0, 2], [2, 1]])
    ends = np.array([[2, 3], [3, 2]])
    path = ["A", "A", "GD", "A"]

    expected_paths = np.array([["3G", "1G1"], ["1G1", "1G1"]], dtype=object)
    expected_likelihoods = np.array([[0.6, 0.9], [0.9, 0.9]])

    paths, likelihoods = calculate_paths_and_likelihoods(
        branch_lengths, CNs, data, starts, ends, path
    )

    assert np.array_equal(paths, expected_paths)
    assert np.array_equal(likelihoods, expected_likelihoods)


#test_calculate_paths_and_likelihoods()


def add_BP_likelihoods(all_structures, p_up, p_down):

    file = PRECOMPUTED_FILE_FOLDER + "/subbed_mat_u"+str(int(p_up))+"_d"+str(int(p_down))+"_p8_v4.precomputed_paths.pickle"

    data = pkl.load(open(file,'rb'))


    for chrom in all_structures.keys():
        timing_struct_to_BP_likelihood_per_chrom(
            data=data,
            structures=all_structures[chrom],
            p_up=p_up,
            p_down=p_down
            )
    return all_structures

def test_get_BP_likelihoods():
    trees_and_timings = {
        "chrom1": [
            (
                None,
                None,
                None,
                np.array([[2, 3], [3, 2]]),
                None,
            ),
            (
                None,
                None,
                None,
                np.array([[0, 2], [2, 1]]),
                None,
            ),
        ]
    }
    pre, mid, post = 2, 1, 1
    p_up, p_down = 65, 10
    data = {
        "3G": [None, [0.2, 0.6, 0.2]],
        "1G1": [None, [0.9, 0.1, 0]],
        "2G0": [None, [1, 0, 0]],
    }

    expected_BP_likelihoods = {
        "chrom1": [
            np.array([-0.916290731874155, -0.916290731874155]),
            np.array([-0.40546510810816444, -0.40546510810816444]),
        ]
    }

    BP_likelihoods = get_BP_likelihoods(trees_and_timings, pre, mid, post, p_up, p_down, data)

    for chrom in BP_likelihoods.keys():
        for i in range(len(BP_likelihoods[chrom])):
            assert np.allclose(BP_likelihoods[chrom][i], expected_BP_likelihoods[chrom][i], rtol=1e-9, atol=1e-9)



#test_get_BP_likelihoods()



##### STEP 8; from the branch lengths and the BP likelihoods calculate the join CN-SNV likelihoods
##### 
##### 
##### 
##### 
##### 


import numpy as np

def get_poisson_loglikelihood(counts, stacked_branch_lengths, plambda, chrom, lengths, unique_CNs):

    stacked_branch_lengths_float = stacked_branch_lengths.astype(float)

    A = np.log(stacked_branch_lengths_float * plambda * lengths[chrom]) * counts
    A = np.where((np.isnan(A)) & (stacked_branch_lengths_float == 0), -np.inf, A)

    B = -np.tile([scipy.special.gammaln(x + 1) for x in counts], (stacked_branch_lengths.shape[0], 1))

    C = -stacked_branch_lengths_float * plambda * lengths[chrom]

    summed = A + B + C
    total = np.sum(summed, axis=1)

    summed_numeric = summed.astype(float)

    if np.isnan(summed_numeric.flatten()).any():
        logging.warning("#####################\n"*10)
        logging.warning("counts: %s", counts)
        logging.warning("plambda: %s", plambda)
        logging.warning("chrom: %s", chrom)
        logging.warning("lengths: %s", lengths)

        for i in range(stacked_branch_lengths.shape[0]):
            logging.warning("Row %s:", i+1)
            logging.warning("stacked_branch_lengths: %s", stacked_branch_lengths[i])
            logging.warning("A: %s", A[i])
            logging.warning("B: %s", B[i])
            logging.warning("C: %s", C[i])
            logging.warning("summed: %s", summed[i])

        logging.warning("total: %s", total)
        sys.exit()

    # FIX 1, make sure that were not testing if 0 SNVs occured in 0 time:
    # Find the indices where branch_lengths is zero
    zero_indices_branch_lengths = np.argwhere(stacked_branch_lengths == 0)

    # Make sure that the position in summed_numeric corresponding to a zero in branch_lengths and counts is also zero
    for idx in zero_indices_branch_lengths:
        if counts[idx[1]] == 0:
            summed_numeric[idx[0], idx[1]] = 0  # set the specific index to 0


    # FIX 2, make sure that we are not testing to see SNVs on CN 0 chromosomes:
    # If CNs is zero in any column, set every entry in summed_numeric in that column to zero
    for j, cn in enumerate(unique_CNs):
        if cn == 0:
            summed_numeric[:, j] = 0  # set the entire column to 0


    return summed_numeric


def test_get_poisson_loglikelihood():
    counts = np.array([2, 3, 5, 7, 11])
    stacked_branch_lengths = np.array(
        [
            [0.5, 0.6, 0.7, 0.8, 0.9],
            [0.1, 0.2, 0.3, 0.4, 0.5],
        ]
    )
    plambda = 0.2
    chrom = "chrom1"
    to_delete = [1, 3]
    lengths = {"chrom1": 1000}
    expected_total = [np.array([-351.94186892 -125.86288737]), np.array([-477.80475629])]
    total = get_poisson_loglikelihood(counts, stacked_branch_lengths, plambda, chrom, to_delete, lengths)
    assert np.allclose(total, expected_total, rtol=1e-9, atol=1e-9), f"get_poisson_loglikelihood test failed, expected: {expected_total}, but got: {total}"

#test_get_poisson_loglikelihood()

def get_all_poisson_loglikelihoods_per_chr(chrom_structures, plambda, observed_SNV_multiplicities, chrom):
    SNV_likelihoods = []

    for i in range(len(chrom_structures)):
        # get the counts of SNVs per unique CN
        counts = []
        unique_CNs = chrom_structures[i]["unique_CNs"]
        for CN in unique_CNs:
            if int(CN) not in observed_SNV_multiplicities:
                counts += [0]
            else:
                counts += [observed_SNV_multiplicities[CN]]

        assert(unique_CNs == sorted(unique_CNs,reverse=True))

        this_SNV_likelihoods = get_poisson_loglikelihood(
            counts=counts,
            stacked_branch_lengths=chrom_structures[i]["stacked_branch_lengths"],
            plambda=plambda,
            chrom=chrom,
            lengths=LENGTHS,
            unique_CNs=unique_CNs
        )
        
        this_SNV_likelihood = np.sum(this_SNV_likelihoods, axis=1)
        chrom_structures[i]["SNV_log_likelihoods"] = this_SNV_likelihoods #np.transpose(this_SNV_likelihoods)
        #sys.exit()

        chrom_structures[i]["counts"] = counts
        chrom_structures[i]["SNV_log_likelihood"] = this_SNV_likelihood
        chrom_structures[i]["SNV_BP_log_likelihood"] = chrom_structures[i]["SNV_log_likelihood"] + chrom_structures[i]["BP_loglikelihood"]
        SNV_likelihoods += [chrom_structures[i]["SNV_BP_log_likelihood"]]
        chrom_structures[i]["observed_SNV_multiplicities"]  = observed_SNV_multiplicities
    return SNV_likelihoods


def test_get_all_poisson_loglikelihoods_per_chr():
    timings = [
        (None, None, None, [0, 1, 2], None),
        (None, None, None, [0, 2, 4], None)
    ]
    plambda = 0.2
    BP_likelihoods = [np.array([-1, -2, -3]), np.array([-4, -5, -6])]
    observed_SNV_multiplicities = {1: 3, 2: 5}
    chrom = "chrom1"

    expected_SNV_likelihoods = [
        np.array([-17.39516599, -19.73780159, -22.08043718]),
        np.array([-21.13183242, -23.47446761, -25.8171028])
    ]

    SNV_likelihoods = get_all_poisson_loglikelihoods_per_chr(timings, plambda, BP_likelihoods, observed_SNV_multiplicities, chrom)

    for i in range(len(SNV_likelihoods)):
        assert np.allclose(SNV_likelihoods[i], expected_SNV_likelihoods[i], rtol=1e-9, atol=1e-9)



# test_get_all_poisson_loglikelihoods_per_chr()


def find_best_SNV_likelihood(plambda, these_structures, observed_SNV_multiplicities):
    SNV_likelihoods = {}
    best = {}
    chroms = these_structures.keys()

    for chrom in chroms:
        # ERROR SNV_likelihoods[chrom] is saved in these_structures so doesn't need to get passed back like this
        SNV_likelihoods[chrom] = get_all_poisson_loglikelihoods_per_chr(
            chrom_structures = these_structures[chrom],
            plambda = plambda,
            observed_SNV_multiplicities = observed_SNV_multiplicities[chrom],
            chrom=chrom
        )
        
        best[chrom] = (-np.Inf, 0, 0)

        for tree_index in range(len(SNV_likelihoods[chrom])):
            max_SNV_loglik = max(SNV_likelihoods[chrom][tree_index])
            best_timings_row_index = np.argmax(SNV_likelihoods[chrom][tree_index])

            if max_SNV_loglik > best[chrom][0]:
                best[chrom] = (max_SNV_loglik, tree_index, best_timings_row_index)

    total = sum([best[chrom][0] for chrom in chroms])
    return total, best


def test_find_best_SNV_likelihood():
    plambda = 0.2
    timings = {
        "chrom1": [
            (None, None, None, [0, 1, 2], None),
            (None, None, None, [0, 2, 4], None)
        ]
    }
    BP_likelihoods = {
        "chrom1": [np.array([-1, -2, -3]), np.array([-4, -5, -6])]
    }

    expected_total = -17.39516599
    expected_best = {
        "chrom1": (-17.39516599, 0, 0)
    }

    total, best = find_best_SNV_likelihood(plambda, timings, BP_likelihoods,observed_SNV_multiplicities)

    assert np.isclose(total, expected_total, rtol=1e-9, atol=1e-9)

    for chrom in best:
        assert np.isclose(best[chrom][0], expected_best[chrom][0], rtol=1e-9, atol=1e-9)
        assert best[chrom][1:] == expected_best[chrom][1:]

def get_best_struct(best_structure_indicies,best_structures):
    best_structure = {}
    for chrom in best_structures:
        best_tree = best_structure_indicies[chrom][1]
        best_structure[chrom] = copy.deepcopy(best_structures[chrom][best_tree])

        best_row = best_structure_indicies[chrom][2]
        for key in ['epochs_created', 'branch_lengths', 'stacked_branch_lengths', 'starts', 'ends', 'paths', 'BP_individual_log_likelihoods', 'BP_loglikelihood', 'SNV_log_likelihood', 'SNV_BP_log_likelihood', 'SNV_log_likelihoods']:
            best_structure[chrom][key] = best_structure[chrom][key][best_row]

    return(best_structure)


def are_all_chromosomes_viable_by_BP_likelihood(all_structures, tree_flexibility):
    """
    Check if all chromosomes have a viable tree.

    A tree is considered viable if its maximum depth is less than or equal to
    the calculated epochs.

    Args:
    all_structures (dict): A dictionary where each key is a chromosome and its
        value is a list of results. Each result is another dictionary that has
        'tree', 'pre', 'mid' and 'post' keys.

    Returns:
    bool: True if all chromosomes have at least one viable tree, False otherwise.
    """
    return_val = True

    for chrom in all_structures:
        tests = [tree_in_bounds(result["tree"], calculate_epochs(result["pre"],result["mid"],result["post"]) ,tree_flexibility) for result in all_structures[chrom]]
        #tests = [max_tree_depth(result["tree"]) <= calculate_max_allowed_BP_depth(result["pre"], result["mid"], result["post"]) for result in all_structures[chrom]]
        if not any(tests):
            return_val = False

    if not return_val:
        logging.info("print what failed")
        logging.info("all_structures %s", all_structures)
        result = all_structures[0][0]
        max_epoch = calculate_epochs(result["pre"],result["mid"],result["post"])
        logging.info("max_epochs %s", max_epoch)
        for chrom in all_structures:
            logging.info("chrom %s", chrom)
            for result in all_structures[chrom]: 
                logging.info("print the trees: %s", result["tree"] )
                logging.info("prihnt the tree depths %s", max_tree_depth(result["tree"]))
            tests = [tree_in_bounds(result["tree"], calculate_epochs(result["pre"],result["mid"],result["post"]) ,tree_flexibility) for result in all_structures[chrom]]
            logging.info("tests %s", tests)
            logging.info("any tests %s", any(tests))
            logging.info("tree_flexibility")
            logging.info(f"DEPTH >= {max_epoch + 3 - tree_flexibility} and depth <= {max_epoch + 3}")


    return return_val



def find_BP_and_SNV_loglik(plambda_start, p_up_start, p_down_start, p_window, plambda_window, all_structures, observed_SNV_multiplicities, total_SNVs, tree_flexibility):
    best_neg_loglik = float("inf")
    best_p_up = 0
    best_p_down = 0
    best_plambda = 0
    best_structures = None
    for p_up in range(max(0,p_up_start-p_window), min(100,p_up_start+p_window+1)): #[p_up_start]:
        for p_down in range(max(0,p_down_start-p_window), min(100,p_down_start+p_window+1)): 
            if p_up+p_down>100:
                continue
            # get the branching process loglik here, 
            # then can pass it in to optimise over the SNVs, otherwise it gets recomputed so many times
            add_BP_likelihoods(all_structures, p_up, p_down)
            temp_structures = copy.deepcopy(all_structures)
            
            def optimize_func(plambda):
                if total_SNVs > 0:
                    assert plambda != 0, f"plambda: {plambda}, observed_SNV_multiplicities: {observed_SNV_multiplicities}"

                total, best = find_best_SNV_likelihood(
                    plambda=plambda,
                    these_structures=all_structures,
                    observed_SNV_multiplicities=observed_SNV_multiplicities
                    )
                return -total

            bounds = (plambda_start * plambda_window, plambda_start / plambda_window)
            x0 = plambda_start
            res = opt.minimize_scalar(optimize_func, bounds=bounds, method='bounded', options={'disp': True})
            #for chrom in all_structures:
            #    for result in all_structures[chrom]:
            #        has_non_inf_values = np.logical_not(np.isinf(result["SNV_BP_log_likelihood"])).any()
            #        if not has_non_inf_values:
            #            print("chrom:", chrom)
            #            print(result)

            if res.fun < best_neg_loglik:
                best_neg_loglik = res.fun
                best_p_up = p_up
                best_p_down = p_down
                best_plambda = res.x
                best_structures = temp_structures #copy.deepcopy(all_structures)
                best_res = res

    # need to work out why this is getting None when it shouldn't
    if best_structures == None:
        # find out why is it none:
        # is it because that for one chromosome each of the tree structures are too deep? check it:
        if not are_all_chromosomes_viable_by_BP_likelihood(all_structures, tree_flexibility):
            print("ERROR, all chroms have viable trees, need to implement code and investigate why the SNVs are throwing this off.")
            sys.exit()
        return None
    total, best_structure_indicies = find_best_SNV_likelihood(
        plambda=best_plambda,
        these_structures=best_structures,
        observed_SNV_multiplicities=observed_SNV_multiplicities
    )
    assert best_neg_loglik == -total
    best_structure = get_best_struct(best_structure_indicies,best_structures)

    return best_neg_loglik, best_p_up, best_p_down, best_plambda, best_structure


def test_find_BP_and_SNV_loglik():
    plambda_start = 0.2
    p_up_start = 50
    p_down_start = 50
    trees_and_timings = {
        "chrom1": [
            (None, None, None, [0, 1, 2], None),
            (None, None, None, [0, 2, 4], None)
        ]
    }
    pre = 0
    mid = 0
    post = 0
    p_window = 0
    plambda_window = 1

    expected_loglik = 17.39516599
    expected_best_p_up = 50
    expected_best_p_down = 50
    expected_best_plambda = 0.2

    loglik, best_p_up, best_p_down, best_plambda, res = find_BP_and_SNV_loglik(plambda_start, p_up_start, p_down_start, trees_and_timings, pre, mid, post, p_window, plambda_window, observed_SNV_multiplicities)

    assert np.isclose(loglik, expected_loglik, rtol=1e-9, atol=1e-9)
    assert best_p_up == expected_best_p_up
    assert best_p_down == expected_best_p_down
    assert np.isclose(best_plambda, expected_best_plambda, rtol=1e-9, atol=1e-9)



#test_find_BP_and_SNV_loglik()


##### 
##### 
##### 

def sum_SNV_counts(observed_SNV_multiplicities):
    d = observed_SNV_multiplicities
    total = 0
    for key1 in d:
        for key2 in d[key1]:
            total += d[key1][key2]
    return total


def sum_chrom_multiplicities(observed_CN_multiplicities):
    return sum(observed_CN_multiplicities.values())


def test_sum_SNV_counts():
    # Test case 1
    observed_SNV_multiplicities1 = {
        "chr1": {0: 5, 1: 10},
        "chr2": {1: 7, 2: 3},
    }
    expected_total1 = 25
    assert sum_SNV_counts(observed_SNV_multiplicities1) == expected_total1

    # Test case 2
    observed_SNV_multiplicities2 = {
        "chr1": {0: 3, 1: 4},
        "chr2": {2: 1, 3: 2},
    }
    expected_total2 = 10
    assert sum_SNV_counts(observed_SNV_multiplicities2) == expected_total2



test_sum_SNV_counts()


def test_sum_chrom_multiplicities():
    # Test case 1
    observed_CN_multiplicities1 = {"chr1": 5, "chr2": 10}
    expected_total1 = 15
    assert sum_chrom_multiplicities(observed_CN_multiplicities1) == expected_total1

    # Test case 2
    observed_CN_multiplicities2 = {"chr1": 3, "chr2": 7}
    expected_total2 = 10
    assert sum_chrom_multiplicities(observed_CN_multiplicities2) == expected_total2



test_sum_chrom_multiplicities()


test_sum_chrom_multiplicities()


def sum_observed_CNs(dictionary):
    total_sum = 0
    for value_list in dictionary.values():
        total_sum += sum(value_list)
    return total_sum

def test_sum_observed_CNs():
    test_cases = [
        {
            "input": {0: [1, 0], 1: [2, 0], 2: [2, 2]},
            "expected_output": 7,
        },
        {
            "input": {3: [3, 0], 4: [1, 0], 5: [1, 0]},
            "expected_output": 5,
        },
        {
            "input": {6: [1, 0], 7: [5, 1], 8: [3, 3]},
            "expected_output": 13,
        },
        {
            "input": {9: [2, 1], 10: [1, 0], 11: [1, 0]},
            "expected_output": 5,
        },
    ]

    for test_case in test_cases:
        input_dictionary = test_case["input"]
        expected_output = test_case["expected_output"]
        result = sum_observed_CNs(input_dictionary)
        assert result == expected_output, f"For {input_dictionary}, expected {expected_output} but got {result}"


test_sum_observed_CNs()




def print_results(res, path, p_up, p_down, pre, mid, post):
    logging.info("############################\n" * 10)
    logging.info(res)
    logging.info("path: " + str(path))
    logging.info("p_up: " + str(p_up))
    logging.info("p_down: " + str(p_down))
    logging.info("pre: " + str(pre))
    logging.info("mid: " + str(mid))
    logging.info("post: " + str(post))


def print_all_trees_and_timings(trees_and_timings):
    logging.debug("Print all trees and timings info")
    for chrom in trees_and_timings:
        logging.debug(chrom)
        logging.debug(trees_and_timings[chrom])

def initial_rate_estimate(pre, mid, post, SS):
    logging.info("Initial rate estimate")
    total_SNV_time = pre + mid + post + 2 + 1
    total_SNVs = sum_SNV_counts(SS["observed_SNV_multiplicities"])
    logging.info("total SNVs")
    logging.info(total_SNVs)
    logging.info("observed_SNV_multiplicities")
    logging.info(SS["observed_SNV_multiplicities"])

    logging.info("observed_CNs")
    logging.info(SS["observed_CNs"])

    total_chromosomes = sum_observed_CNs(SS["observed_CNs"])
    plambda_start = float(total_SNVs) / float(total_SNV_time) / float(total_chromosomes) * 23
    logging.info("plambda_start: " + str(plambda_start))

    return plambda_start, total_SNVs

def get_tree_and_rate_parameters(res, SEARCH_DEPTH, SS):
    if res == SEARCH_DEPTH:
        p_up = int(SS["p_up"] * 100)
        p_down = int(SS["p_down"] * 100)
        pre = SS["pre"]
        mid = SS["mid"]
        post = SS["post"]
        path = generate_path(pre, mid, post)
    else:
        path = SS["searchable_likelihoods"]["path"].iloc[res]
        p_up = int(SS["searchable_likelihoods"]['p_up'].iloc[res])
        p_down = int(SS["searchable_likelihoods"]['p_down'].iloc[res])
        pre, mid, post = path_code_to_pre_mid_post(path)



    for x in [pre, mid, post, p_up, p_down]:
        assert isinstance(x, (int, np.integer)), f"Expected integer, but got {x} of type {type(x)}."

    for x in [p_up,p_down]:
        assert(x >= 0 and x <= 100)

    return path, p_up, p_down, pre, mid, post

import math

def poisson_zero_probability(lam,rounds):
    """
    Calculate the probability that a Poisson random variable is zero for a given parameter.

    :param lam: The Poisson parameter (average rate).
    :return: The probability that a Poisson random variable is zero.
    """
    num_SNVs = 0
    probability = (math.exp(-lam) * (lam ** num_SNVs)) / math.factorial(num_SNVs)
    probability = (probability**23)**rounds
    return probability

from scipy.stats import dirichlet

def calculate_log_likelihood(alpha, probabilities):
    """
    Calculate the log-likelihood of given probabilities based on the Dirichlet distribution.

    Args:
        alpha: list or np.array
            Parameters of the Dirichlet distribution.
        probabilities: np.array
            Probabilities to calculate the log-likelihood for. 

    Returns:
        float: The log-likelihood of the probabilities based on the Dirichlet distribution.
    """
    # Reverse the scaling and rounding transformation
    probabilities = [float(x)/100 for x in probabilities]
    print(probabilities)

    return dirichlet.logpdf(probabilities, alpha)

def test_calculate_log_likelihood():
    alpha = [20, 20, 100]
    probabilities = generate_dirichlet_probability(alpha)*100
    log_likelihood = calculate_log_likelihood(alpha, probabilities)

    assert isinstance(log_likelihood, np.float64), f"Output is not a float: {type(log_likelihood)}"

test_calculate_log_likelihood()


def handle_results(result_dict):
    for key, value in result_dict.items():
        logging.info(f'{key}: {value}')
        if isinstance(key,int):
            handle_chroms(result_dict, key)

def handle_chroms(result_dict, chrom):
    for x in ["observed_SNV_multiplicities","tree","CNs","starts","ends","branch_lengths","paths","BP_individual_log_likelihoods","BP_loglikelihood","unique_CNs","stacked_branch_lengths","counts"]:
        logging.info(f'{x}: {result_dict[chrom][x]}')
    temp = [count / LENGTHS[chrom] for count in result_dict[chrom]["counts"]]
    #print(temp)
    #print(result_dict[chrom]["stacked_branch_lengths"])
    #print(result_dict[chrom])
    normalize_counts(result_dict, chrom, temp)
    handle_likelihoods(result_dict, chrom)

def normalize_counts(result_dict, chrom, temp):
    for x in range(len(temp)):
        if result_dict[chrom]["stacked_branch_lengths"][x] != 0:
            temp[x] = temp[x]/result_dict[chrom]["stacked_branch_lengths"][x]
    logging.info(f'normalised_counts: {temp}')

def handle_likelihoods(result_dict, chrom):
    for x in ["SNV_log_likelihoods","SNV_log_likelihood","SNV_BP_log_likelihood"]:
        logging.info(f'{x}: {result_dict[chrom][x]}')

def check_failure(all_structures):
    SNV_failure = False
    for chrom in all_structures:
        for structure in all_structures[chrom]:
            if np.isinf(structure["SNV_log_likelihood"]).all():
                print("chrom",chrom)
                print(all_structures[chrom])
                print("failed because of SNV likelihood")
                SNV_failure=True
    return SNV_failure

import time



def handle_errors(all_structures):
    print("failed and don't know why")
    SNV_failure = check_failure(all_structures)
    if SNV_failure:
        print("isn't a SNV failure")
        return True

    for i in range(20):
        print("SEEMINGLY EXCLUSIVELY A SNV OR A BP FAILURE")
    for chrom in all_structures:
        print("chrom", chrom)
        for structure in all_structures[chrom]:
            for key in structure:
                print(f"\tKEY:{key}, {structure[key]}")
    for i in range(20):
        print("***************************************")

    impossible = False
    for chrom in all_structures:
        this_chrom_possible = False
        for structure in all_structures[chrom]:
            if any(str(x) != '-inf' for x in structure["SNV_BP_log_likelihood"]):
                this_chrom_possible = True

        if not this_chrom_possible:
            impossible = True

    return impossible


def find_solutions(SS, p_window, plambda_window, tree_flexibility, alpha):
    SS["observed_SNV_multiplicities"] = count_SNV_multiplicities(SS['simulated_chromosomes'])
    SEARCH_DEPTH = len(SS['searchable_likelihoods'])
    results = []

    aggregated_execution_times = {
        "total_time": 0,
        "get_all_trees_and_timings": 0,
        "timing_struct_to_all_structures": 0,
        "find_BP_and_SNV_loglik": 0
    }

    for res in range(SEARCH_DEPTH + 1):
        start_time_res = time.time()
        print("RES", res)
        path, p_up_start, p_down_start, pre, mid, post = get_tree_and_rate_parameters(res, SEARCH_DEPTH, SS)
        print_results(res, path, p_up_start, p_down_start, pre, mid, post)

        total_epochs = calculate_epochs(pre, mid, post)
        assert total_epochs == pre + mid + post + 2
        max_epoch = total_epochs

        start_time_trees_and_timings = time.time()
        trees = get_all_trees(
            observed_SNV_multiplicities=SS["observed_SNV_multiplicities"],
            observed_CNs=SS["observed_CNs"],
            total_epochs=total_epochs,
            tree_flexibility=tree_flexibility
        )

        if trees is not None:
            trees_and_timings = get_all_timings(
                trees=trees,
                observed_SNV_multiplicities=SS["observed_SNV_multiplicities"],
                total_epochs=total_epochs,
                pre=pre,
                mid=mid,
                post=post
            )
        else:
            trees_and_timings = None

        end_time_trees_and_timings = time.time()

        aggregated_execution_times["get_all_trees_and_timings"] += round(end_time_trees_and_timings - start_time_trees_and_timings, 2)

        if trees_and_timings is None:
            continue

        start_time_all_structures = time.time()
        all_structures = timing_struct_to_all_structures(
            trees_and_timings=trees_and_timings,
            pre=pre,
            mid=mid,
            post=post,
            max_epoch=max_epoch
        )
        end_time_all_structures = time.time()

        aggregated_execution_times["timing_struct_to_all_structures"] += round(end_time_all_structures - start_time_all_structures, 2)

        plambda_start, total_SNVs = initial_rate_estimate(pre, mid, post, SS)

        start_time_BP_SNV_loglik = time.time()
        BP_SNV_loglik = find_BP_and_SNV_loglik(
            plambda_start=plambda_start,
            p_up_start=p_up_start,
            p_down_start=p_down_start,
            p_window=p_window,
            plambda_window=plambda_window,
            all_structures=all_structures,
            observed_SNV_multiplicities=SS["observed_SNV_multiplicities"],
            total_SNVs=total_SNVs,
            tree_flexibility=tree_flexibility
        )
        end_time_BP_SNV_loglik = time.time()

        aggregated_execution_times["find_BP_and_SNV_loglik"] += round(end_time_BP_SNV_loglik - start_time_BP_SNV_loglik, 2)

        if BP_SNV_loglik is not None:
            best_neg_loglik, best_p_up, best_p_down, best_plambda, best_structure = BP_SNV_loglik
            logging.info("BP_SNV_output")

            probabilities = [best_p_up, best_p_down, 100-best_p_up-best_p_down]
            log_likelihood = calculate_log_likelihood(alpha, probabilities)

            result_dict = copy.deepcopy(best_structure)
            result_dict["est_neg_loglik"] = best_neg_loglik-log_likelihood
            result_dict["est_p_up"] = best_p_up
            result_dict["est_p_down"] = best_p_down
            result_dict["est_plambda"] = best_plambda
        else:
            result_dict = None

        if result_dict is not None:
            handle_results(result_dict)
        else:
            impossible = handle_errors(all_structures)
            if impossible:
                continue

        result_dict["execution_times"] = aggregated_execution_times
        results.append(result_dict)

    return results




def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Run all processing scripts with the given test case and optional parameters."
    )
    
    parser.add_argument(
        "-t", "--test_case", 
        type=int, 
        default=0, 
        help="The name or number of the test case you want to process."
    )
    
    parser.add_argument(
        "-b", "--debug", 
        action='store_true', 
        help="Enable debug logging."
    )
    
    parser.add_argument(
        "-f", "--simulation_filename", 
        type=str, 
        default="simulation", 
        help="The name of the simulation and analysis files."
    )
    
    parser.add_argument(
        "-p", "--p_window", 
        type=int, 
        default=1, 
        help="Specify the p_window value."
    )
    
    parser.add_argument(
        "-l", "--plambda_window", 
        type=float, 
        default=0.01, 
        help="Specify the plambda_window value."
    )
    
    parser.add_argument(
        "-x", "--tree_flexibility", 
        type=int, 
        default=4, 
        help="Specify the tree flexibility value."
    )
    
    parser.add_argument(
        "-a", "--alpha", 
        nargs='+', 
        type=float, 
        default=[2, 2, 10],
        help="Alpha parameter for Dirichlet distribution (optional, default is [2, 2, 10])."
    )

    return parser.parse_args()


def init(args):
    test_case = args.test_case
    SS = load_results_from_file(test_case=test_case, simulation_filename=args.simulation_filename)
    return SS



def save_dict_to_file(dictionary, test_case, simulation_filename):
    filename = f'{SIMULATIONS_FILE_FOLDER}/{simulation_filename}_{test_case}.pickle'
    
    if os.path.exists(filename):
        with open(filename, 'rb') as f:
            old_data = pickle.load(f)
        old_data.update(dictionary)
        dictionary = old_data
    
    with open(filename, 'wb') as f:
        pickle.dump(dictionary, f)



def main():
    args = parse_arguments()
    print(args)

    if args.debug:
        logging.basicConfig(level = logging.DEBUG)
    else:
        logging.basicConfig(level = logging.INFO)

    logging.info("load and start finding results")

    SS = init(args)

    SS["solutions"] = find_solutions(SS=SS, p_window=args.p_window, plambda_window=args.plambda_window, tree_flexibility=args.tree_flexibility, alpha=args.alpha)

    # Sort the dictionary based on 'best_neg_loglik' in descending order
    sorted_solutions = sorted(SS['solutions'], key=lambda x: x['est_neg_loglik'], reverse=True)

    print(f"the truth is:")
    print(
        f"pre: {SS['pre']},\n"
        f"mid: {SS['mid']},\n"
        f"post: {SS['post']},\n"
        f"p_up: {SS['p_up']},\n"
        f"p_down: {SS['p_down']},\n"
        f"rate: {SS['rate']}"
    )
    for solution in sorted_solutions:
        print(
            f"best_neg_loglik: {solution['est_neg_loglik']}, "
            f"best_p_up: {solution['est_p_up']}, "
            f"best_p_down: {solution['est_p_down']}, "
            f"best_plambda: {solution['est_plambda']}, "
            f"pre: {solution[0]['pre']}, "
            f"mid: {solution[0]['mid']}, "
            f"post: {solution[0]['post']}, \n"
            f"total_time: {solution['execution_times']['total_time']}, "
            f"time_get_all_trees_and_timings: {solution['execution_times']['get_all_trees_and_timings']}, "
            f"time_timing_struct_to_all_structures: {solution['execution_times']['timing_struct_to_all_structures']}, "
            f"time_find_BP_and_SNV_loglik: {solution['execution_times']['timing_struct_to_all_structures']} \n\n"
        )

   
    for chrom in SS["CN_trees"]:
        print(chrom, SS["CN_trees"][chrom])
    print(SS.keys())

    # Save some solutions to file here
    save_dict_to_file(dictionary=SS, test_case=args.test_case, simulation_filename=args.simulation_filename)

def load_results_from_file(test_case, simulation_filename):
    file_name = f'{SIMULATIONS_FILE_FOLDER}/{simulation_filename}_{test_case}.pickle'
    with open(file_name, 'rb') as f:
        data = pickle.load(f)
    return data


if __name__ == "__main__":
    main()

