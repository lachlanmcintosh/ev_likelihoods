
def test_insert_node_into_tree_structure():
    # Test case 1
    tree1 = {
        "unique_identifier": 1,
        "child": None,
        "complement": None,
        "epoch_created": 1
    }
    node1 = {
        "unique_identifier": 2,
        "parent": 1,
        "epoch_created": 2
    }
    expected_tree1 = {
        "unique_identifier": 1,
        "child": {
            "unique_identifier": 2,
            "parent": 1,
            "epoch_created": 2,
            "child": None,
            "complement": None
        },
        "complement": {
            "unique_identifier": 1,
            "child": None,
            "complement": None,
            "epoch_created": 2
        },
        "epoch_created": 1
    }
    result_tree1 = insert_node_into_tree_structure(tree1, node1)
    assert result_tree1 == expected_tree1, f"Test case 1 failed: Expected:\n{json.dumps(expected_tree1, indent=4)}\nBut got:\n{json.dumps(result_tree1, indent=4)}"

    # Test case 2
    tree2 = {
        "unique_identifier": 1,
        "child": {
            "unique_identifier": 2,
            "parent": 1,
            "child": None,
            "complement": None,
            "epoch_created": 2
        },
        "complement": None,
        "epoch_created": 1
    }
    node2 = {
        "unique_identifier": 3,
        "parent": 2,
        "epoch_created": 3
    }
    expected_tree2 = {
        "unique_identifier": 1,
        "child": {
            "unique_identifier": 2,
            "parent": 1,
            "child": {
                "unique_identifier": 3,
                "parent": 2,
                "epoch_created": 3,
                "child": None,
                "complement": None
            },
            "complement": {
                "unique_identifier": 2,
                "parent": 1,
                "child": None,
                "complement": None,
                "epoch_created": 3
            },
            "epoch_created": 2
        },
        "complement": None,
        "epoch_created": 1
    }
    result_tree2 = insert_node_into_tree_structure(tree2, node2)
    assert result_tree2 == expected_tree2, f"Test case 2 failed: Expected:\n{json.dumps(expected_tree2, indent=4)}\nBut got:\n{json.dumps(result_tree2, indent=4)}"

    # Test case 3
    tree3 = None
    node3 = {
        "unique_identifier": 1,
        "parent": None,
        "epoch_created": 2
    }
    expected_tree3 = None
    result_tree3 = insert_node_into_tree_structure(tree3, node3)
    assert result_tree3 == expected_tree3, f"Test case 3 failed: Expected:\n{json.dumps(expected_tree3, indent=4)}\nBut got:\n{json.dumps(result_tree3, indent=4)}"

    # Test case 4
    tree4 = {
        "unique_identifier": 1,
        "child": None,
        "complement": None
    }
    node4 = None
    expected_tree4 = {
        "unique_identifier": 1,
        "child": None,
        "complement": None
    }
    result_tree4 = insert_node_into_tree_structure(tree4, node4)
    assert result_tree4 == expected_tree4, f"Test case 4 failed: Expected:\n{json.dumps(expected_tree4, indent=4)}\nBut got:\n{json.dumps(result_tree4, indent=4)}"


test_insert_node_into_leaf()
test_insert_node_under_complement()
test_insert_node_into_tree_structure()



def test_add_copynumber_to_tree_structure():
    # Test case 1: Simple tree with no children or complements
    tree = {
        "child": None,
        "complement": None,
        "dead": False,
    }
    expected_tree = {
        "child": None,
        "complement": None,
        "dead": False,
        "copy_number": 1,
    }
    result_tree = add_copynumber_to_tree_structure(tree)
    assert result_tree == expected_tree, f"Expected {expected_tree}, but got {result_tree}"

    # Test case 2: Tree with one level of children
    tree = {
        "child": {
            "child": None,
            "complement": None,
            "dead": False,
        },
        "complement": {
            "child": None,
            "complement": None,
            "dead": True,
        },
    }
    expected_tree = {
        "child": {
            "child": None,
            "complement": None,
            "dead": False,
            "copy_number": 1,
        },
        "complement": {
            "child": None,
            "complement": None,
            "dead": True,
            "copy_number": 0,
        },
        "copy_number": 1,
    }
    result_tree = add_copynumber_to_tree_structure(tree)
    assert result_tree == expected_tree, f"Expected {expected_tree}, but got {result_tree}"

    # Test case 3: Tree with multiple levels
    tree = {
        "child": {
            "child": {
                "child": None,
                "complement": None,
                "dead": False,
            },
            "complement": {
                "child": None,
                "complement": None,
                "dead": True,
            },
        },
        "complement": {
            "child": None,
            "complement": None,
            "dead": False,
        },
    }
    expected_tree = {
        "child": {
            "child": {
                "child": None,
                "complement": None,
                "dead": False,
                "copy_number": 1,
            },
            "complement": {
                "child": None,
                "complement": None,
                "dead": True,
                "copy_number": 0,
            },
            "copy_number": 1,
        },
        "complement": {
            "child": None,
            "complement": None,
            "dead": False,
            "copy_number": 1,
        },
        "copy_number": 2,
    }
    result_tree = add_copynumber_to_tree_structure(tree)
    assert result_tree == expected_tree, f"Expected {expected_tree}, but got {result_tree}"

test_add_copynumber_to_tree_structure()

