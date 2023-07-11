from typing import Any, Dict, List, Optional, Tuple, Union

def convert_dict_tree_to_list(
    tree: Optional[Dict[str, Any]],
    total_epochs: Optional[int] = None,
    is_truth: bool = False
) -> Optional[List[Union[Tuple[int, int], List]]]:
    """
    Convert a dictionary tree to a list.

    Args:
        tree: The dictionary tree to convert.
        total_epochs: The total number of epochs if given.
        is_truth: Flag to indicate if the truth value is considered.

    Returns:
        A list representation of the tree, or None if the input tree is None.
    """
    if tree is None:
        return None

    copy_number = tree.get('copy_number')
    epoch_created = tree.get('epoch_created')

    if is_truth:
        epoch_created = _get_epoch_created(tree, total_epochs)

    child_tree = convert_dict_tree_to_list(tree.get('child'), total_epochs, is_truth)
    complement_tree = convert_dict_tree_to_list(tree.get('complement'), total_epochs, is_truth)

    if child_tree is None:
        return [(copy_number, epoch_created)]
    return [(copy_number, epoch_created), child_tree, complement_tree]

def _get_epoch_created(tree: Dict[str, Any], total_epochs: Optional[int]) -> int:
    """
    Determines the epoch_created based on the presence of 'child' key in the tree.

    Args:
        tree: The dictionary tree to check.
        total_epochs: The total number of epochs.

    Returns:
        An integer value representing the epoch created.
    """
    if 'child' in tree and tree["child"] is not None:
        return tree["child"]["epoch_created"]
    else:
        assert total_epochs is not None, "total_epochs is not provided but required when 'child' key is not present."
        return total_epochs

def get_dicts_list(
    trees_list: List[Any],
    total_epochs: int,
    is_truth: bool
) -> List[Dict[str, Any]]:
    """
    Get a list of dictionaries representing trees.

    Args:
        trees_list: List of trees to convert.
        total_epochs: Total number of epochs.
        is_truth: Whether the truth value is considered.

    Returns:
        A list of dictionaries representing trees.
    """
    dicts_list = []
    for tree in trees_list:
        list_representation = convert_dict_tree_to_list(tree, total_epochs, is_truth)
        dicts_list.append(list_representation)
    return dicts_list

