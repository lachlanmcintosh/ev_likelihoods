import copy
import pytest
from insert_node import insert_node_into_leaf, insert_node_under_complement

@pytest.mark.parametrize("tree,node,expected", [
    (
        {'unique_identifier': 1, 'child': None, 'complement': None},
        {'unique_identifier': 2, 'parent': 1, 'epoch_created': 100},
        {'unique_identifier': 1,
        'child': {'unique_identifier': 2, 'parent': 1, 'epoch_created': 100, 'child': None, 'complement': None},
        'complement': {'unique_identifier': 1, 'child': None, 'complement': None, 'epoch_created': 100}}
    ),
    # ... add more test cases as needed ...
])
def test_insert_node_into_leaf(tree, node, expected):
    insert_node_into_leaf(tree, node)
    assert tree == expected, f'Expected {expected}, but got {tree}'

@pytest.mark.parametrize("tree,node,expected_error", [
    (
        {'unique_identifier': 1, 'child': {}, 'complement': None},
        {'unique_identifier': 2, 'parent': 1, 'epoch_created': 100},
        "Expected tree['child'] to be None"
    ),
    (
        {'unique_identifier': 1, 'child': None, 'complement': {}},
        {'unique_identifier': 2, 'parent': 1, 'epoch_created': 100},
        "Expected tree['complement'] to be None"
    ),
    # ... add more test cases as needed ...
])
def test_insert_node_into_leaf_errors(tree, node, expected_error):
    with pytest.raises(AssertionError, match=expected_error):
        insert_node_into_leaf(tree, node)

# for insert_node_under_complement

@pytest.mark.parametrize("tree,node,expected", [
    # typical case
    (
        {'unique_identifier': 1, 'child': {'unique_identifier': 2, 'parent': 1, 'child': None, 'complement': None}, 'complement': {'unique_identifier': 1, 'child': None, 'complement': None}},
        {'unique_identifier': 3, 'parent': 1, 'epoch_created': 1, 'SNVs': [], 'paternal': False, 'dead': False},
        {'unique_identifier': 1,
         'child': {'unique_identifier': 2, 'parent': 1, 'child': None, 'complement': None},
         'complement': {'unique_identifier': 1, 
                        'child': {'unique_identifier': 3, 'parent': 1, 'epoch_created': 1, 'child': None, 'complement': None, 'SNVs': [], 'paternal': False, 'dead': False},
                        'complement': {'unique_identifier': 1, 'child': None, 'complement': None}}
        }
    ),
    # ... add more test cases as needed ...
])
def test_insert_node_under_complement(tree, node, expected):
    insert_node_under_complement(tree, node)
    assert tree == expected, f"""Expected:\n{json.dumps(expected, indent=4)}\nBut got:\n{json.dumps(tree, indent=4)}"""

@pytest.mark.parametrize("tree,node", [
    (
        {'unique_identifier': 1, 'child': {'unique_identifier': 2, 'parent': 1, 'child': None, 'complement': None}, 'complement': {'unique_identifier': 1, 'child': None, 'complement': None}},
        {'unique_identifier': 3, 'parent': 1, 'epoch_created': 1, 'extra_key': 'extra_value', 'SNVs': [], 'paternal': False, 'dead': False},
    ),
    # ... add more test cases as needed ...
])
def test_insert_node_under_complement_errors(tree, node):
    with pytest.raises(AssertionError):
        insert_node_under_complement(tree, node)

# for insert_node_into_tree_structure

@pytest.mark.parametrize("tree,node,expected", [
    (
        {'unique_identifier': 1, 'child': None, 'complement': None, 'epoch_created': 1},
        {'unique_identifier': 2, 'parent': 1, 'epoch_created': 2},
        {'unique_identifier': 1, 
         'child': {'unique_identifier': 2, 'parent': 1, 'epoch_created': 2, 'child': None, 'complement': None}, 
         'complement': {'unique_identifier': 1, 'child': None, 'complement': None, 'epoch_created': 2}, 
         'epoch_created': 1}
    ),
    # ... add more test cases as needed ...
])
def test_insert_node_into_tree_structure(tree, node, expected):
    result_tree = insert_node_into_tree_structure(tree, node)
    assert result_tree == expected, f"""Expected:\n{json.dumps(expected, indent=4)}\nBut got:\n{json.dumps(result_tree, indent=4)}"""

@pytest.mark.parametrize("tree,node", [
    (
        {'unique_identifier': 1, 'child': None, 'complement': None, 'epoch_created': 1},
        {'unique_identifier': None, 'parent': 1, 'epoch_created': 2},
    ),
    (
        {'unique_identifier': None, 'child': None, 'complement': None, 'epoch_created': 1},
        {'unique_identifier': 2, 'parent': 1, 'epoch_created': 2},
    ),
])
def test_insert_node_into_tree_structure_errors(tree, node):
    with pytest.raises(KeyError):
        insert_node_into_tree_structure(tree, node)
