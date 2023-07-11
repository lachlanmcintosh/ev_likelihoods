


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


