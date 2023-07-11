import pytest
from remove_redundancy import remove_SNVs_from_tree_structure, remove_dead_leaf_nodes_from_tree_structure, remove_redundant_parents_from_tree_structure

# Tree fixture
@pytest.fixture(params=[
    {'unique_identifier': 0, 'parent': -1, 'epoch_created': 0, 'paternal': True, 'dead': False, 'child': None, 'complement': None, 'copy_number': 1, 'SNV_multiplicity': 0},
    {'dead': True, 'parent': 1, 'child': None, 'complement': None},
    {'parent': 1, 'copy_number': 1, 'child': {'parent': 1, 'copy_number': 1, 'child': None, 'complement': None, 'epoch_created': 2}, 'complement': {'parent': 1, 'copy_number': 0, 'child': None, 'complement': None, 'epoch_created': 2}, 'epoch_created': 1}
])
def tree(request):
    return request.param

# Parametrized test cases for remove_SNVs_from_tree_structure
@pytest.mark.parametrize("tree,expected", [
    (tree[0], tree[0].copy()),
    (tree[1], None)
])
def test_remove_SNVs_from_tree_structure(tree, expected):
    result = remove_SNVs_from_tree_structure(tree)
    assert result == expected, f"remove_SNVs_from_tree_structure({tree}) = {result}, expected {expected}"

# Parametrized test cases for remove_dead_leaf_nodes_from_tree_structure
@pytest.mark.parametrize("tree,expected", [
    (None, None),
    (tree[0], tree[0].copy()),
    (tree[1], None)
])
def test_remove_dead_leaf_nodes_from_tree_structure(tree, expected):
    result = remove_dead_leaf_nodes_from_tree_structure(tree)
    assert result == expected, f"remove_dead_leaf_nodes_from_tree_structure({tree}) = {result}, expected {expected}"

# Parametrized test cases for remove_redundant_parents_from_tree_structure
@pytest.mark.parametrize("node,expected", [
    (None, None),
    (tree[0], tree[0].copy()),
    (tree[2], {'parent': 1, 'copy_number': 1, 'child': None, 'complement': None, 'epoch_created': 1})
])
def test_remove_redundant_parents_from_tree_structure(node, expected):
    result = remove_redundant_parents_from_tree_structure(node)
    assert result == expected, f"remove_redundant_parents_from_tree_structure({node}) = {result}, expected {expected}"
