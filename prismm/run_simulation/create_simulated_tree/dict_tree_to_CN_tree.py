from typing import Dict, List, Tuple

from clonal_trees.utils.make_left_heavy import make_left_heavy

def get_child_and_complement_trees_from_truth_tree(truth_tree: Dict) -> Tuple:
    """
    Get child and complement trees from a given truth tree.

    :param truth_tree: The input truth tree.
    :return: The child and complement trees.
    """
    child_tree = convert_truth_tree_to_CN_tree(truth_tree["child"]) if truth_tree.get("child") else None
    complement_tree = convert_truth_tree_to_CN_tree(truth_tree["complement"]) if truth_tree.get("complement") else None
    return child_tree, complement_tree

def convert_truth_tree_to_CN_tree(truth_tree: Dict) -> List:
    """
    Convert truth tree to copy number (CN) tree.

    :param truth_tree: The input truth tree.
    :return: The resultant CN tree.
    """
    child_tree, complement_tree = get_child_and_complement_trees_from_truth_tree(truth_tree)
    
    CN_tree = [truth_tree.get("copy_number", 0)]
    if child_tree:
        CN_tree.append(child_tree)
    if complement_tree:
        CN_tree.append(complement_tree)

    return CN_tree

def convert_truth_trees_to_CN_trees(truth_trees: Dict[str, Dict]) -> Dict[str, List]:
    """
    Convert truth trees to copy number (CN) trees and make the trees left-heavy.

    :param truth_trees: The input truth trees.
    :return: The resultant CN trees.
    """
    for chrom_type in truth_trees:
        truth_trees[chrom_type] = convert_truth_tree_to_CN_tree(truth_trees[chrom_type])
        truth_trees[chrom_type] = make_left_heavy(truth_trees[chrom_type])
    return truth_trees
