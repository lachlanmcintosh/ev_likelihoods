
def test_add_SNV_multiplicity_to_tree_structure():
    # Test case 1: simple tree with no children or complements
    tree = {
        "copy_number": 1,
        "child": None,
        "complement": None,
        "epoch_created": 0,
        "SNVs": [{"epoch_created": 0}, {"epoch_created": 1}]
    }
    expected_tree = tree.copy()
    expected_tree["SNV_multiplicity"] = 2
    result_tree = add_SNV_multiplicity_to_tree_structure(tree)
    assert result_tree == expected_tree, f"Expected {expected_tree}, but got {result_tree}"

    # Test case 2: tree with one level of children
    tree = {
        "copy_number": 2,
        "epoch_created": 0,
        "SNVs": [{"epoch_created": 0}],
        "child": {
            "copy_number": 1,
            "epoch_created": 1,
            "SNVs": [{"epoch_created": 1}],
            "child": None,
            "complement": None
        },
        "complement": {
            "copy_number": 1,
            "epoch_created": 1,
            "SNVs": [{"epoch_created": 2}],
            "child": None,
            "complement": None
        }
    }
    expected_tree = {
        "copy_number": 2,
        "epoch_created": 0,
        "SNV_multiplicity": 1,
        "SNVs": [{"epoch_created": 0}],
        "child": {
            "copy_number": 1,
            "epoch_created": 1,
            "SNV_multiplicity": 1,
            "SNVs": [{"epoch_created": 1}],
            "child": None,
            "complement": None
        },
        "complement": {
            "copy_number": 1,
            "epoch_created": 1,
            "SNV_multiplicity": 1,
            "SNVs": [{"epoch_created": 2}],
            "child": None,
            "complement": None
        }
    }
    result_tree = add_SNV_multiplicity_to_tree_structure(tree)
    assert result_tree == expected_tree, f"Expected {expected_tree}, but got {result_tree}"

    # Additional test cases can be added for more complex trees

test_add_SNV_multiplicity_to_tree_structure()


def test_remove_SNVs_from_tree_structure():
    # Test case 1
    tree = {
        "SNVs": ["dummy"],
        "child": None,
        "complement": None
    }
    expected_tree = {
        "child": None,
        "complement": None
    }
    result_tree = remove_SNVs_from_tree_structure(tree)
    assert result_tree == expected_tree, f"Expected {expected_tree}, but got {result_tree}"



def test_remove_dead_leaf_nodes_from_tree_structure():
    # Test case 1
    tree = {
        "dead": True,
        "parent": 1,
        "child": None,
        "complement": None
    }
    expected_tree = None
    result_tree = remove_dead_leaf_nodes_from_tree_structure(tree)
    assert result_tree == expected_tree, f"Expected {expected_tree}, but got {result_tree}"



def test_are_all_descendants_zero_in_tree_structure():
    # Test case 1
    node = {
        "copy_number": 0,
        "child": None,
        "complement": None
    }
    expected_result = True
    result = are_all_descendants_zero_in_tree_structure(node)
    assert result == expected_result, f"Expected {expected_result}, but got {result}"

    # Add more test cases as needed

test_remove_SNVs_from_tree_structure()
test_remove_dead_leaf_nodes_from_tree_structure()
test_are_all_descendants_zero_in_tree_structure()


def test_remove_redundant_parents_from_tree_structure():
    # Test case 1
    node = {
        "parent": 1,
        "copy_number": 1,
        "child": {
            "parent": 1,
            "copy_number": 1,
            "child": None,
            "complement": None,
            "epoch_created": 2
        },
        "complement": {
            "parent": 1,
            "copy_number": 0,
            "child": None,
            "complement": None,
            "epoch_created": 2
        },
        "epoch_created": 1
    }
    expected_node = {
        "parent": 1,
        "copy_number": 1,
        "child": None,
        "complement": None,
        "epoch_created": 1
    }
    result_node = remove_redundant_parents_from_tree_structure(node)
    assert result_node == expected_node, f"Expected {expected_node}, but got {result_node}"

    # Add more test cases as needed

test_remove_redundant_parents_from_tree_structure()



def test_assign_epoch_killed_to_tree_structure():
    # Test case 1
    tree = {
        "epoch_created": 1,
        "child": {
            "epoch_created": 2,
            "child": None,
            "complement": None
        },
        "complement": {
            "epoch_created": 2,
            "child": None,
            "complement": None
        }
    }
    max_epochs = 3
    expected_tree = {
        "epoch_created": 1,
        "epoch_killed": 2,
        "child": {
            "epoch_created": 2,
            "epoch_killed": 3,
            "child": None,
            "complement": None
        },
        "complement": {
            "epoch_created": 2,
            "epoch_killed": 3,
            "child": None,
            "complement": None
        }
    }
    result_tree = assign_epoch_killed_to_tree_structure(tree, max_epochs)
    assert result_tree == expected_tree, f"Expected {expected_tree}, but got {result_tree}"

    # Add more test cases as needed

test_assign_epoch_killed_to_tree_structure()



def test_convert_truth_tree_to_CN_tree():
    # Case 1: Test with a simple truth_tree
    truth_tree = {
        "copy_number": 1,
        "child": {
            "copy_number": 2,
            "child": None,
            "complement": None
        },
        "complement": {
            "copy_number": 3,
            "child": None,
            "complement": None
        }
    }
    expected_CN_tree = [1, [2], [3]]
    result_CN_tree = convert_truth_tree_to_CN_tree(truth_tree)
    assert result_CN_tree == expected_CN_tree, f"Expected {expected_CN_tree}, but got {result_CN_tree}"

test_convert_truth_tree_to_CN_tree()



def test_convert_truth_trees_to_CN_trees():
    # Case 1: Test with simple truth_trees
    truth_trees = {
        "type1": {
            "copy_number": 1,
            "child": {
                "copy_number": 2,
                "child": None,
                "complement": None
            },
            "complement": {
                "copy_number": 3,
                "child": None,
                "complement": None
            }
        }
    }
    expected_CN_trees = {
        "type1": [1, [3], [2]]
    }
    result_CN_trees = convert_truth_trees_to_CN_trees(truth_trees)
    assert result_CN_trees == expected_CN_trees, f"Expected {expected_CN_trees}, but got {result_CN_trees}"

test_convert_truth_trees_to_CN_trees()


def test_count_copy_numbers():
    simulated_chromosomes = {
        "type1": [
            {"paternal": True, "dead": False},
            {"paternal": True, "dead": True},
            {"paternal": False, "dead": False},
            {"paternal": False, "dead": False},
        ],
        "type2": [
            {"paternal": True, "dead": False},
            {"paternal": True, "dead": False},
            {"paternal": False, "dead": True},
            {"paternal": False, "dead": True},
        ],
    }

    expected_counts = {
        "type1": [1, 2],  # 1 living paternal and 2 living non-paternal for type1
        "type2": [2, 0],  # 2 living paternal and 0 living non-paternal for type2
    }

    observed_counts = count_copy_numbers(simulated_chromosomes)

    assert observed_counts == expected_counts, f"Expected {expected_counts}, but got {observed_counts}"

test_count_copy_numbers()


def test_count_CN_multiplicities():
    observed_CNs = {
        "type1": [1, 2, 2, 3, 3, 3],
        "type2": [2, 2, 3, 3, 4],
    }

    expected_multiplicities = {
        1: 1,
        2: 4,
        3: 5,
        4: 1,
    }

    observed_multiplicities = count_CN_multiplicities(observed_CNs)

    assert observed_multiplicities == expected_multiplicities, f"Expected {expected_multiplicities}, but got {observed_multiplicities}"


def test_count_copy_number_multiplicities():
    observed_copy_numbers = {
        "type1": [1, 2, 2, 3, 3, 3],
        "type2": [2, 2, 3, 3, 4],
    }

    expected_multiplicities = {
        1: 1,
        2: 4,
        3: 5,
        4: 1,
    }

    observed_multiplicities = count_copy_number_multiplicities(observed_copy_numbers)

    assert observed_multiplicities == expected_multiplicities, f"Expected {expected_multiplicities}, but got {observed_multiplicities}"


test_count_CN_multiplicities()
test_count_copy_number_multiplicities()


def test_generate_dirichlet_probability():
    alpha = [20.0, 20.0, 100.0] # [0.5, 0.5, 0.5]
    probabilities = generate_dirichlet_probability(alpha)

    assert isinstance(probabilities, np.ndarray), f"Output is not a numpy array: {type(probabilities)}"
    assert np.isclose(probabilities.sum(), 1, atol=1e-2), f"Probabilities do not sum to 1 within tolerance: {probabilities.sum()}"
    assert (probabilities >= 0).all() and (probabilities <= 1).all(), f"Probabilities contain values outside [0, 1]: {probabilities}"
    assert np.isclose(np.round(probabilities, 2).sum(), 1, atol=1e-8), "Rounded probabilities do not sum to 1"
    #assert np.round(probabilities, 2).sum() == 1, "Rounded probabilities do not sum to 1"


test_generate_dirichlet_probability()

from scipy.stats import poisson
from dataclasses import dataclass


