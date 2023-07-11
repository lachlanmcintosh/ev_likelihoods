import copy
from typing import Dict, List, Tuple, Optional

from clonal_trees.run_simulation.create_simulated_tree.dict_tree_to_CN_tree import convert_truth_tree_to_CN_tree, get_child_and_complement_trees_from_truth_tree
from clonal_trees.utils.make_left_heavy import make_left_heavy


def insert_node_into_leaf(tree: Dict, node: Dict) -> None:
    """
    Inserts a node into a leaf of a given tree.

    :param tree: The tree where the node is to be inserted.
    :param node: The node to be inserted into the tree.
    """

    # Assert that the given tree is a leaf node
    assert tree["child"] is None, f"Expected tree['child'] to be None, but got {tree['child']}"
    assert tree["complement"] is None, f"Expected tree['complement'] to be None, but got {tree['complement']}"

    # Create a deep copy of the tree as the complement
    tree["complement"] = copy.deepcopy(tree)
    tree["complement"]["epoch_created"] = node["epoch_created"]

    # Insert the node as a child
    tree["child"] = copy.deepcopy(node)
    tree["child"]["child"] = None
    tree["child"]["complement"] = None



def insert_node_under_complement(tree: Dict, node: Dict) -> None:
    """
    Inserts a node under the complement of a given tree.

    :param tree: The tree where the node is to be inserted.
    :param node: The node to be inserted under the complement of the tree.
    """

    # nodes are inserted into the tree in their order of their unique identifier
    # which maps to the order of chromosomes created.
    # there may be many children of "complement" which is why it gets dragged down the tree
    # to find out when it bifurcated into each new chromosome it is the parent of

    # complement cannot be None, because the "child" node is always created first
    # and complement is always created with it as a place holder
    assert tree["complement"] is not None, f"Expected tree['complement'] to not be None, but it was."

    if node["parent"] == -1:
        # these can't get passed further on because they are the original chromosomes
        expected_keys = {"unique_identifier", "epoch_created", "parent", "SNVs", "paternal", "dead"}
        assert set(node.keys()) == expected_keys, f"Unexpected keys: {node.keys()} - Expected keys: {expected_keys}"

        assert node["unique_identifier"] < 46, f"Expected node['unique_identifier'] to be less than 46, but got {node['unique_identifier']}"

        tree["complement"] = copy.deepcopy(node)
        tree["complement"]["child"] = None
        tree["complement"]["complement"] = None
    else:
        # we don't know how far down the complement line that we have to insert the node: so just keep going until it is done:

        # insert_node_under_complement and insert_node_under_leaf do not return anything,
        # they modify the tree structure given to them,
        # why do we return something here from insert_node_into_tree_structure? CHECKHERE

        tree["complement"] = insert_node_into_tree_structure(tree["complement"], node)

import copy
from typing import Dict, Optional

def insert_node_into_tree(child_tree: Optional[Dict], complement_tree: Optional[Dict], node: Optional[Dict]) -> Tuple[Optional[Dict], Optional[Dict]]:
    """
    Insert a node into the child and complement trees if they are not None.

    :param child_tree: The child tree where the node is to be inserted.
    :param complement_tree: The complement tree where the node is to be inserted.
    :param node: The node to be inserted into the trees.
    :return: The modified child and complement trees.
    """
    if child_tree is not None:
        child_tree = insert_node_into_tree_structure(child_tree, node)
    if complement_tree is not None:
        complement_tree = insert_node_into_tree_structure(complement_tree, node)

    return child_tree, complement_tree

def insert_node_into_tree_structure(tree: Optional[Dict], node: Optional[Dict]) -> Optional[Dict]:
    """
    Inserts a node into a tree structure according to the unique identifier and parent-child relationships.

    :param tree: The tree where the node is to be inserted.
    :param node: The node to be inserted into the tree.
    :return: The modified tree.
    """
    if node is None or tree is None:
        return tree

    if node["unique_identifier"] is not None and tree["unique_identifier"] is not None and node["unique_identifier"] != tree["unique_identifier"]:
        if node["parent"] == tree["unique_identifier"]:
            if tree["child"] is None:
                insert_node_into_leaf(tree, node)
            else:
                insert_node_under_complement(tree, node)
        else:
            tree["child"], tree["complement"] = insert_node_into_tree(child_tree=tree["child"], complement_tree=tree["complement"], node=node)

    return tree
