from typing import Dict, Optional

def add_copynumber_to_tree_structure(tree: Optional[Dict]) -> Optional[Dict]:
    """
    Add the "copy_number" attribute to each node in the given tree structure. 
    This function works by recursively adding up the copy number of the leaves up to the root.

    :param tree: The tree to which the "copy_number" attribute is to be added.
    :return: The modified tree with added "copy_number" attribute.
    """
    if "copy_number" in tree:
        return tree

    if tree["child"] is None and tree["complement"] is None:
        tree["copy_number"] = 0 if tree.get("dead") else 1
    else:
        if tree["child"] is not None:
            tree["child"] = add_copynumber_to_tree_structure(tree["child"])
        if tree["complement"] is not None:
            tree["complement"] = add_copynumber_to_tree_structure(tree["complement"])

        child_copy_number = tree["child"]["copy_number"] if tree["child"] else 0
        complement_copy_number = tree["complement"]["copy_number"] if tree["complement"] else 0
        tree["copy_number"] = child_copy_number + complement_copy_number

    return tree

def add_SNV_multiplicity_to_tree_structure(tree: Optional[Dict]) -> Optional[Dict]:
    """
    Adds "SNV_multiplicity" attribute to each node in the given tree structure.

    :param tree: The tree to which the "SNV_multiplicity" attribute is to be added.
    :return: The modified tree with added "SNV_multiplicity" attribute.
    """
    if tree["copy_number"] == 0:
        tree["SNV_multiplicity"] = None
        return tree

    if tree["child"] is None:
        assert tree["complement"] is None
        count = sum(1 for SNV in tree["SNVs"] if SNV["epoch_created"] >= tree["epoch_created"])
        tree["SNV_multiplicity"] = count
        return tree

    tree["child"] = add_SNV_multiplicity_to_tree_structure(tree["child"])
    tree["complement"] = add_SNV_multiplicity_to_tree_structure(tree["complement"])

    count = sum(1 for SNV in tree["SNVs"] if tree["epoch_created"] <= SNV["epoch_created"] < tree["child"]["epoch_created"])
    tree["SNV_multiplicity"] = count

    return tree

def assign_epoch_killed_to_tree_structure(tree: Optional[Dict], max_epochs: int) -> Optional[Dict]:
    """
    The function calculates the epoch_killed for each node in the given tree.

    :param tree: A dictionary representing a tree node with keys:
                 'child', 'complement', 'epoch_created', and 'epoch_killed'.
    :param max_epochs: An integer representing the maximum epoch that epoch_killed can be.
    :return: The modified tree with 'epoch_killed' values added.
    """
    if tree is None:
        return None

    child = tree.get('child')
    complement = tree.get('complement')

    if child is not None or complement is not None:
        child_epoch_created = child.get('epoch_created')
        complement_epoch_created = complement.get('epoch_created')

        if child_epoch_created != complement_epoch_created:
            raise ValueError("Epoch created values of child and complement do not match.")

        tree['child'] = assign_epoch_killed_to_tree_structure(child, max_epochs)
        tree['complement'] = assign_epoch_killed_to_tree_structure(complement, max_epochs)

        tree['epoch_killed'] = child['epoch_created']
    else:
        tree['epoch_killed'] = max_epochs
        
    return tree
