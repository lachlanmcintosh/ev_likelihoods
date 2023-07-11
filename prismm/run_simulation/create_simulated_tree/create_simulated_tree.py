import logging
from typing import Dict

from clonal_trees.run_simulation.create_simulated_tree.insert_node import insert_node_into_tree_structure
from clonal_trees.run_simulation.create_simulated_tree.add_metadata import (
    add_copynumber_to_tree_structure, 
    add_SNV_multiplicity_to_tree_structure, 
    assign_epoch_killed_to_tree_structure
)
from clonal_trees.run_simulation.create_simulated_tree.remove_redundancy import (
    remove_SNVs_from_tree_structure, 
    remove_dead_leaf_nodes_from_tree_structure, 
    remove_redundant_parents_from_tree_structure
)

def create_root_node() -> Dict:
    """
    Create the root node of the tree.

    :return: The root node of the tree.
    """
    return {
        'unique_identifier': -1,
        'parent': None,
        'epoch_created': 0,  # TODO; ERROR HERE; this should be None since it isn't a real node
        'paternal': None,
        'child': None,
        'complement': None,
        'SNVs': []
    }

def insert_nodes_into_tree(tree: Dict, sorted_list: list) -> Dict:
    """
    Insert all nodes and add metadata to tree and simplify.

    :param tree: The tree into which nodes are to be inserted.
    :param sorted_list: The sorted list of nodes to be inserted.
    :return: The tree with nodes inserted.
    """
    for new_node in sorted_list:
        tree = insert_node_into_tree_structure(tree, new_node[1])  
        # new_node[1] is the node itself, new_node[0] is the unique_identifier  
    return tree

def create_truth_trees_from_simulation(simulated_chromosomes: Dict, max_epochs: int) -> Dict:
    """
    Given the simulated chromosomes, the function creates the truth trees.
    It creates a separate tree for each type of chromosome.
    Each tree is updated and simplified by adding metadata to the nodes,
    adding copy numbers, SNV multiplicity, removing SNVs and dead nodes,
    removing redundant nodes, and assigning 'epoch_killed' to nodes.

    :param simulated_chromosomes: A dictionary of simulated chromosomes.
    :param max_epochs: The maximum number of epochs for simulation.
    :return: A dictionary of truth trees for each chromosome type.
    """
    trees = {}
    min_identifier = min(node['unique_identifier'] for chrom_type in simulated_chromosomes for node in simulated_chromosomes[chrom_type])
    assert min_identifier == 0, "Expected the minimum 'unique_identifier' to be 0."

    for chrom_type in simulated_chromosomes:
        sorted_list = sorted([(x["unique_identifier"],x) for x in simulated_chromosomes[chrom_type]])
        tree = create_root_node()
        tree = insert_nodes_into_tree(tree, sorted_list)

        tree = add_copynumber_to_tree_structure(tree)
        tree = add_SNV_multiplicity_to_tree_structure(tree)
        tree = remove_SNVs_from_tree_structure(tree)
        tree = remove_dead_leaf_nodes_from_tree_structure(tree)
        tree = remove_redundant_parents_from_tree_structure(tree)
        tree = assign_epoch_killed_to_tree_structure(tree, max_epochs)

        trees[chrom_type] = tree

    return trees
