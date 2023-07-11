from typing import Optional, Dict

def are_all_descendants_zero_in_tree_structure(node: Optional[Dict]) -> bool:
    """
    Checks if the "copy_number" attribute of the given node and all its descendants in the tree structure is zero.

    :param node: The node to start the check from.
    :return: True if all "copy_number" attributes are zero, False otherwise.
    """
    if node is None:
        return True

    copy_number = node.get("copy_number")
    if copy_number is None or copy_number != 0:
        return False

    return are_all_descendants_zero_in_tree_structure(node.get("child")) and are_all_descendants_zero_in_tree_structure(node.get("complement"))
        
def check_and_remove_redundant_node_in_tree_structure(node: Optional[Dict], main_child: Optional[Dict], other_child: Optional[Dict]) -> Optional[Dict]:
    """
    Checks if a node in the tree structure is redundant based on its children's "copy_number" attributes,
    and if it is, it updates the "epoch_created" of the main child and removes the node.

    :param node: The node to check.
    :param main_child: The main child node of the node.
    :param other_child: The other child node of the node.
    :return: The updated node (if it was redundant, it gets replaced with the main child node).
    """
    if main_child is not None:
        main_child_copy_number = main_child.get("copy_number")
        node_copy_number = node.get("copy_number")

        if main_child_copy_number is not None and node_copy_number is not None and main_child_copy_number == node_copy_number:
            if not are_all_descendants_zero_in_tree_structure(other_child):
                raise ValueError("Invalid tree: child and complement copy_number do not add up to parent's copy_number")

            # Update epoch_created of the main_child to that of its parent
            main_child["epoch_created"] = node.get("epoch_created")
            return main_child

    return node
