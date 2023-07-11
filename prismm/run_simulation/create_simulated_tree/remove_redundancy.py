from typing import Optional, Dict

from clonal_trees.run_simulation.create_simulated_tree.check_trees import check_and_remove_redundant_node_in_tree_structure

def remove_SNVs_from_tree_structure(tree: Optional[Dict]) -> Optional[Dict]:
    """
    Recursively removes the "SNVs" attribute from each node in the given tree structure.

    :param tree: The tree from which the "SNVs" attribute is to be removed.
    :return: The modified tree without the "SNVs" attribute.
    """
    if tree is None:
        return None

    tree.pop("SNVs", None)

    if tree.get("child"):
        tree["child"] = remove_SNVs_from_tree_structure(tree["child"])
    if tree.get("complement"):
        tree["complement"] = remove_SNVs_from_tree_structure(tree["complement"])

    return tree

def remove_dead_leaf_nodes_from_tree_structure(tree: Optional[Dict]) -> Optional[Dict]:
    """
    Recursively removes dead leaf nodes from the given tree structure.

    :param tree: The tree from which dead leaf nodes are to be removed.
    :return: The modified tree with dead leaf nodes removed.
    """
    if tree is None:
        return None

    tree['child'] = remove_dead_leaf_nodes_from_tree_structure(tree['child'])
    tree['complement'] = remove_dead_leaf_nodes_from_tree_structure(tree['complement'])

    if tree['child'] is None and tree['complement'] is None and tree.get('dead', False) and tree['parent'] != -1:
        return None

    return tree
    
def remove_redundant_parents_from_tree_structure(node: Optional[Dict]) -> Optional[Dict]:
    """
    Recursively removes redundant parents from a tree structure. A node is considered redundant if its "copy_number"
    equals its child's "copy_number" and all the descendants of the other child node have a "copy_number" of zero.

    :param node: The root of the tree.
    :return: The root of the tree after redundant nodes have been removed.
    """
    if node is None:
        return None

    # Recursively process child and complement nodes
    if node.get("child"):
        node["child"] = remove_redundant_parents_from_tree_structure(node["child"])
    if node.get("complement"):
        node["complement"] = remove_redundant_parents_from_tree_structure(node["complement"])

    # Check for the conditions to remove the node, only if the node has a parent
    if node["parent"] is not None:
        node = check_and_remove_redundant_node_in_tree_structure(node, node.get("child"), node.get("complement"))
        node = check_and_remove_redundant_node_in_tree_structure(node, node.get("complement"), node.get("child"))

    return node
