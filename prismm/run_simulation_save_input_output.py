
from constants import LENGTHS
import numpy as np
import random
import copy
import sys
import shelve

def count_paternity(chromosomes, paternal):
    with open('IO_LOGGER/count_paternity_io_log.txt', 'a') as f:
        f.write('Input: {0}\n'.format(repr((chromosomes, paternal))))
    "\n    Count the number of chromosomes of a particular chromosome type with a specific paternal type that are not marked as dead.\n\n    Args:\n        chromosomes (list): A list of dictionaries representing chromosomes of a particular chromosomal type/number. Each chromosome dictionary \n            must have a 'paternal' key indicating the type of paternal chromosome and a 'dead' key \n            indicating whether the chromosome is marked as dead or not.\n        paternal (str): The type of paternal chromosome to count.\n\n    Raises:\n        ValueError: If a chromosome is missing the 'paternal' key or the 'dead' key.\n\n    Returns:\n        int: The number of chromosomes of the specified type with the specified paternal type that are not marked as dead.\n    "
    for chrom in chromosomes:
        if ('paternal' not in chrom):
            raise ValueError("Chromosome is missing 'paternal' key")
        if ('dead' not in chrom):
            raise ValueError("Chromosome is missing 'dead' key")
    _result = len([chrom for chrom in chromosomes if ((chrom['paternal'] == paternal) and (not chrom['dead']))])
    with open('IO_LOGGER/count_paternity_io_log.txt', 'a') as f:
        f.write('Output: {0}\n'.format(repr(_result)))
    return len([chrom for chrom in chromosomes if ((chrom['paternal'] == paternal) and (not chrom['dead']))])

pass

def check_all_chrs_are_unique(simulated_chromosomes):
    with open('IO_LOGGER/check_all_chrs_are_unique_io_log.txt', 'a') as f:
        f.write('Input: {0}\n'.format(repr(simulated_chromosomes)))
    '\n    Check if all chromosomes have unique identifiers.\n    '
    ids = [chrom.get('unique_identifier') for chrom_type in simulated_chromosomes.values() for chrom in chrom_type]
    _result = (len(ids) == len(set(ids)))
    with open('IO_LOGGER/check_all_chrs_are_unique_io_log.txt', 'a') as f:
        f.write('Output: {0}\n'.format(repr(_result)))
    return (len(ids) == len(set(ids)))

pass

def check_expected_keys_in_simulated_chromosomes_present(simulated_chromosomes):
    with open('IO_LOGGER/check_expected_keys_in_simulated_chromosomes_present_io_log.txt', 'a') as f:
        f.write('Input: {0}\n'.format(repr(simulated_chromosomes)))
    '\n    Check if all keys that are expected to be present in simulated chromosomes are actually present.\n    '
    expected_keys = ['SNVs', 'paternal', 'epoch_created', 'parent', 'unique_identifier', 'dead']
    for chrom_type in simulated_chromosomes.values():
        for chrom in chrom_type:
            if (sorted(expected_keys) != sorted(chrom.keys())):
                return False
    _result = True
    with open('IO_LOGGER/check_expected_keys_in_simulated_chromosomes_present_io_log.txt', 'a') as f:
        f.write('Output: {0}\n'.format(repr(_result)))
    return True

pass

def get_ev_string(pre, mid, post):
    with open('IO_LOGGER/get_ev_string_io_log.txt', 'a') as f:
        f.write('Input: {0}\n'.format(repr((pre, mid, post))))
    '\n    Generate a list of characters representing the evolutionary path the cancer takes.\n    '
    assert (isinstance(pre, int) and (pre >= (- 1))), 'Invalid pre value'
    assert (isinstance(mid, int) and (mid >= (- 1))), 'Invalid mid value'
    assert (isinstance(post, int) and (post >= (- 1))), 'Invalid post value'
    if (pre == (- 1)):
        assert (mid == (- 1))
        assert (post == (- 1))
    if (mid == (- 1)):
        assert (post == (- 1))
    ev_str = []
    ev_str += (['A'] * pre)
    ev_str += (['G'] * (mid > (- 1)))
    ev_str += (['A'] * mid)
    ev_str += (['G'] * (post > (- 1)))
    ev_str += (['A'] * post)
    _result = ev_str
    with open('IO_LOGGER/get_ev_string_io_log.txt', 'a') as f:
        f.write('Output: {0}\n'.format(repr(_result)))
    return ev_str

pass

def initialize_simulated_chromosomes():
    with open('IO_LOGGER/initialize_simulated_chromosomes_io_log.txt', 'a') as f:
        f.write('Input: {0}\n'.format(repr(())))
    simulated_chromosomes = {}
    for chrom_type in range(23):
        simulated_chromosomes[chrom_type] = [{'unique_identifier': (chrom_type + x), 'parent': (- 1), 'epoch_created': 0, 'paternal': (x == 0), 'SNVs': [], 'dead': False} for x in (0, 23)]
    _result = simulated_chromosomes
    with open('IO_LOGGER/initialize_simulated_chromosomes_io_log.txt', 'a') as f:
        f.write('Output: {0}\n'.format(repr(_result)))
    return simulated_chromosomes

pass

def simulate_snvs(simulated_chromosomes, LENGTHS, rate, epoch):
    with open('IO_LOGGER/simulate_snvs_io_log.txt', 'a') as f:
        f.write('Input: {0}\n'.format(repr((simulated_chromosomes, LENGTHS, rate, epoch))))
    SNV_count = 0
    for chrom_type in simulated_chromosomes:
        for chrom in simulated_chromosomes[chrom_type]:
            if chrom['dead']:
                continue
            additional_SNV_count = np.random.poisson((rate * LENGTHS[chrom_type]), 1)[0]
            for x in range((SNV_count + 1), ((SNV_count + additional_SNV_count) + 1)):
                chrom['SNVs'] += [{'unique_identifier': str(x), 'epoch_created': epoch}]
            SNV_count = (SNV_count + additional_SNV_count)
    _result = (SNV_count, simulated_chromosomes)
    with open('IO_LOGGER/simulate_snvs_io_log.txt', 'a') as f:
        f.write('Output: {0}\n'.format(repr(_result)))
    return (SNV_count, simulated_chromosomes)

pass

def simulate_gd(simulated_chromosomes, epoch, chrom_count):
    with open('IO_LOGGER/simulate_gd_io_log.txt', 'a') as f:
        f.write('Input: {0}\n'.format(repr((simulated_chromosomes, epoch, chrom_count))))
    for chrom_type in simulated_chromosomes:
        new_chromosomes = []
        for chrom in simulated_chromosomes[chrom_type]:
            if chrom['dead']:
                continue
            new_chromosome = copy.deepcopy(chrom)
            chrom_count += 1
            new_chromosome['unique_identifier'] = chrom_count
            new_chromosome['epoch_created'] = epoch
            new_chromosome['parent'] = chrom['unique_identifier']
            new_chromosomes += [new_chromosome]
        simulated_chromosomes[chrom_type] += new_chromosomes
    _result = (chrom_count, simulated_chromosomes)
    with open('IO_LOGGER/simulate_gd_io_log.txt', 'a') as f:
        f.write('Output: {0}\n'.format(repr(_result)))
    return (chrom_count, simulated_chromosomes)

pass

def simulate_anueploidy_agnostic(simulated_chromosomes, epoch, chrom_count):
    with open('IO_LOGGER/simulate_anueploidy_agnostic_io_log.txt', 'a') as f:
        f.write('Input: {0}\n'.format(repr((simulated_chromosomes, epoch, chrom_count))))
    for chrom_type in simulated_chromosomes:
        chroms = simulated_chromosomes[chrom_type]
        for which in random.sample(range(len(chroms)), random.randint(0, len(chroms))):
            chroms[which]['dead'] = True
            chrom_count -= 1
        for which in random.choices(range(len(chroms)), k=random.randint(0, len(chroms))):
            if (not chroms[which]['dead']):
                new_chrom = copy.deepcopy(chroms[which])
                new_chrom['epoch_created'] = epoch
                new_chrom['parent'] = chroms[which]['unique_identifier']
                chrom_count += 1
                new_chrom['unique_identifier'] = chrom_count
                chroms.append(new_chrom)
        simulated_chromosomes[chrom_type] = chroms
    _result = (chrom_count, simulated_chromosomes)
    with open('IO_LOGGER/simulate_anueploidy_agnostic_io_log.txt', 'a') as f:
        f.write('Output: {0}\n'.format(repr(_result)))
    return (chrom_count, simulated_chromosomes)

pass

def simulate_anueploidy(simulated_chromosomes, epoch, chrom_count, p_up, p_down):
    with open('IO_LOGGER/simulate_anueploidy_io_log.txt', 'a') as f:
        f.write('Input: {0}\n'.format(repr((simulated_chromosomes, epoch, chrom_count, p_up, p_down))))
    for chrom_type in simulated_chromosomes:
        while True:
            new_chromosomes = []
            for chrom in simulated_chromosomes[chrom_type]:
                change = np.random.choice([1, 0, (- 1)], 1, p=[p_up, ((1 - p_up) - p_down), p_down])
                if chrom['dead']:
                    new_chromosomes += [chrom]
                    continue
                if (change == 0):
                    new_chromosomes += [chrom]
                elif (change == 1):
                    new_chromosome = copy.deepcopy(chrom)
                    chrom_count += 1
                    new_chromosome['unique_identifier'] = chrom_count
                    new_chromosome['epoch_created'] = epoch
                    new_chromosome['parent'] = chrom['unique_identifier']
                    new_chromosomes += [new_chromosome]
                    new_chromosomes += [chrom]
                elif (change == (- 1)):
                    new_chromosome = copy.deepcopy(chrom)
                    new_chromosome['dead'] = True
                    new_chromosomes += [new_chromosome]
            if (len([x for x in new_chromosomes if (not x['dead'])]) != 0):
                break
        simulated_chromosomes[chrom_type] = new_chromosomes
    _result = (chrom_count, simulated_chromosomes)
    with open('IO_LOGGER/simulate_anueploidy_io_log.txt', 'a') as f:
        f.write('Output: {0}\n'.format(repr(_result)))
    return (chrom_count, simulated_chromosomes)

pass

def check_simulated_chromosomes(simulated_chromosomes, pre, mid, post, ev_sequence):
    with open('IO_LOGGER/check_simulated_chromosomes_io_log.txt', 'a') as f:
        f.write('Input: {0}\n'.format(repr((simulated_chromosomes, pre, mid, post, ev_sequence))))
    assert ((((pre + mid) + post) + 2) == len(ev_sequence))
    if ((post == 0) or ((mid == 0) and (post == (- 1)))):
        assert (ev_sequence[(- 1)] == 'G')
        for chrom_type in simulated_chromosomes:
            if ((count_paternity(simulated_chromosomes[chrom_type], paternal=True) % 2) == 0):
                pretty_print(simulated_chromosomes[chrom_type])
                pretty_print(count_paternity(simulated_chromosomes[chrom_type], paternal=True))
                assert ((count_paternity(simulated_chromosomes[chrom_type], paternal=True) % 2) == 0)
            if ((count_paternity(simulated_chromosomes[chrom_type], paternal=False) % 2) == 0):
                pretty_print(simulated_chromosomes[chrom_type])
                pretty_print(count_paternity(simulated_chromosomes[chrom_type], paternal=False))
                assert ((count_paternity(simulated_chromosomes[chrom_type], paternal=False) % 2) == 0)
    for chrom_type in simulated_chromosomes:
        assert (len(simulated_chromosomes) != 0)
    assert check_all_chrs_are_unique(simulated_chromosomes)
    assert check_expected_keys_in_simulated_chromosomes_present(simulated_chromosomes)
    _result = None
    with open('IO_LOGGER/check_simulated_chromosomes_io_log.txt', 'a') as f:
        f.write('Output: {0}\n'.format(repr(_result)))

def simulate_single_with_poisson_timestamps_names(p_up, p_down, pre, mid, post, rate, agnostic=False):
    with open('IO_LOGGER/simulate_single_with_poisson_timestamps_names_io_log.txt', 'a') as f:
        f.write('Input: {0}\n'.format(repr((p_up, p_down, pre, mid, post, rate, agnostic))))
    SNV_count = 0
    simulated_chromosomes = initialize_simulated_chromosomes()
    chrom_count = 46
    ev_sequence = get_ev_string(pre, mid, post)
    if (len(ev_sequence) == 0):
        return simulated_chromosomes
    for (epoch, epoch_type) in enumerate(ev_sequence):
        (SNV_count, simulated_chromosomes) = simulate_snvs(simulated_chromosomes, LENGTHS, rate, epoch)
        if (((mid != (- 1)) and (epoch == pre)) or ((post != (- 1)) and (epoch == ((pre + 1) + mid)))):
            assert (epoch_type == 'G')
            (chrom_count, simulated_chromosomes) = simulate_gd(simulated_chromosomes, epoch, chrom_count)
        else:
            assert (epoch_type == 'A')
            if agnostic:
                (chrom_count, simulated_chromosomes) = simulate_anueploidy_agnostic(simulated_chromosomes, epoch, chrom_count)
            else:
                (chrom_count, simulated_chromosomes) = simulate_anueploidy(simulated_chromosomes, epoch, chrom_count, p_up, p_down)
    assert ((((pre + mid) + post) + 2) == len(ev_sequence))
    check_simulated_chromosomes(simulated_chromosomes, pre, mid, post, ev_sequence)
    _result = simulated_chromosomes
    with open('IO_LOGGER/simulate_single_with_poisson_timestamps_names_io_log.txt', 'a') as f:
        f.write('Output: {0}\n'.format(repr(_result)))
    return simulated_chromosomes

pass

def insert_node_under_child(tree, node):
    with open('IO_LOGGER/insert_node_under_child_io_log.txt', 'a') as f:
        f.write('Input: {0}\n'.format(repr((tree, node))))
    tree['complement'] = copy.deepcopy(tree)
    tree['complement']['epoch_created'] = node['epoch_created']
    tree['child'] = copy.deepcopy(node)
    tree['child']['child'] = None
    tree['child']['complement'] = None
    _result = None
    with open('IO_LOGGER/insert_node_under_child_io_log.txt', 'a') as f:
        f.write('Output: {0}\n'.format(repr(_result)))

def insert_node_under_complement(tree, node):
    with open('IO_LOGGER/insert_node_under_complement_io_log.txt', 'a') as f:
        f.write('Input: {0}\n'.format(repr((tree, node))))
    if (node['parent'] == (- 1)):
        assert (node['unique_identifier'] < 46)
        tree['complement'] = copy.deepcopy(node)
        tree['complement']['child'] = None
        tree['complement']['complement'] = None
        assert ('copy_number' not in tree['complement'])
    else:
        tree['complement'] = insert_node_into_truth_tree_v2(tree['complement'], node)
    _result = None
    with open('IO_LOGGER/insert_node_under_complement_io_log.txt', 'a') as f:
        f.write('Output: {0}\n'.format(repr(_result)))

def insert_node_into_truth_tree_v2(tree, node):
    with open('IO_LOGGER/insert_node_into_truth_tree_v2_io_log.txt', 'a') as f:
        f.write('Input: {0}\n'.format(repr((tree, node))))
    if ((node is None) or (tree is None)):
        return tree
    if ((node['unique_identifier'] is not None) and (tree['unique_identifier'] is not None) and (node['unique_identifier'] != tree['unique_identifier'])):
        if (node['parent'] == tree['unique_identifier']):
            if (tree['child'] is None):
                insert_node_under_child(tree, node)
            else:
                insert_node_under_complement(tree, node)
        else:
            if (tree['child'] is not None):
                tree['child'] = insert_node_into_truth_tree_v2(tree['child'], node)
            if (tree['complement'] is not None):
                tree['complement'] = insert_node_into_truth_tree_v2(tree['complement'], node)
    _result = tree
    with open('IO_LOGGER/insert_node_into_truth_tree_v2_io_log.txt', 'a') as f:
        f.write('Output: {0}\n'.format(repr(_result)))
    return tree

def add_copynumber_to_tree(tree):
    with open('IO_LOGGER/add_copynumber_to_tree_io_log.txt', 'a') as f:
        f.write('Input: {0}\n'.format(repr(tree)))
    if ('copy_number' in tree):
        return tree
    if ((tree['child'] is None) and (tree['complement'] is None)):
        tree['copy_number'] = (0 if tree.get('dead') else 1)
    else:
        if (tree['child'] is not None):
            tree['child'] = add_copynumber_to_tree(tree['child'])
        if (tree['complement'] is not None):
            tree['complement'] = add_copynumber_to_tree(tree['complement'])
        child_copy_number = (tree['child']['copy_number'] if tree['child'] else 0)
        complement_copy_number = (tree['complement']['copy_number'] if tree['complement'] else 0)
        tree['copy_number'] = (child_copy_number + complement_copy_number)
    _result = tree
    with open('IO_LOGGER/add_copynumber_to_tree_io_log.txt', 'a') as f:
        f.write('Output: {0}\n'.format(repr(_result)))
    return tree

pass

def add_SNV_multiplicity_to_tree(tree):
    with open('IO_LOGGER/add_SNV_multiplicity_to_tree_io_log.txt', 'a') as f:
        f.write('Input: {0}\n'.format(repr(tree)))
    if (tree['copy_number'] == 0):
        tree['SNV_multiplicity'] = None
        return tree
    if (tree['child'] is None):
        assert (tree['complement'] is None)
        count = sum((1 for SNV in tree['SNVs'] if (SNV['epoch_created'] >= tree['epoch_created'])))
        tree['SNV_multiplicity'] = count
        return tree
    assert (tree['child']['epoch_created'] == tree['complement']['epoch_created'])
    tree['child'] = add_SNV_multiplicity_to_tree(tree['child'])
    tree['complement'] = add_SNV_multiplicity_to_tree(tree['complement'])
    count = sum((1 for SNV in tree['SNVs'] if (tree['epoch_created'] <= SNV['epoch_created'] < tree['child']['epoch_created'])))
    tree['SNV_multiplicity'] = count
    _result = tree
    with open('IO_LOGGER/add_SNV_multiplicity_to_tree_io_log.txt', 'a') as f:
        f.write('Output: {0}\n'.format(repr(_result)))
    return tree

pass

def remove_SNVs_from_tree(tree):
    with open('IO_LOGGER/remove_SNVs_from_tree_io_log.txt', 'a') as f:
        f.write('Input: {0}\n'.format(repr(tree)))
    tree.pop('SNVs', None)
    if (tree['child'] is not None):
        assert (tree['complement'] is not None)
        tree['child'] = remove_SNVs_from_tree(tree['child'])
        tree['complement'] = remove_SNVs_from_tree(tree['complement'])
    _result = tree
    with open('IO_LOGGER/remove_SNVs_from_tree_io_log.txt', 'a') as f:
        f.write('Output: {0}\n'.format(repr(_result)))
    return tree

def remove_dead_nodes(tree):
    with open('IO_LOGGER/remove_dead_nodes_io_log.txt', 'a') as f:
        f.write('Input: {0}\n'.format(repr(tree)))
    if (tree is None):
        return None
    tree['child'] = remove_dead_nodes(tree['child'])
    tree['complement'] = remove_dead_nodes(tree['complement'])
    if ((tree['child'] is None) and (tree['complement'] is None) and tree.get('dead', False) and (tree['parent'] != (- 1))):
        return None
    _result = tree
    with open('IO_LOGGER/remove_dead_nodes_io_log.txt', 'a') as f:
        f.write('Output: {0}\n'.format(repr(_result)))
    return tree

def are_all_descendants_zero(node):
    with open('IO_LOGGER/are_all_descendants_zero_io_log.txt', 'a') as f:
        f.write('Input: {0}\n'.format(repr(node)))
    if (node is None):
        return True
    if (node['copy_number'] != 0):
        return False
    _result = (are_all_descendants_zero(node['child']) and are_all_descendants_zero(node['complement']))
    with open('IO_LOGGER/are_all_descendants_zero_io_log.txt', 'a') as f:
        f.write('Output: {0}\n'.format(repr(_result)))
    return (are_all_descendants_zero(node['child']) and are_all_descendants_zero(node['complement']))

pass

pass

pass

def remove_redundant_parents(node):
    with open('IO_LOGGER/remove_redundant_parents_io_log.txt', 'a') as f:
        f.write('Input: {0}\n'.format(repr(node)))

    def check_and_remove_redundant_node(node, main_child, other_child):
        with open('IO_LOGGER/check_and_remove_redundant_node_io_log.txt', 'a') as f:
            f.write('Input: {0}\n'.format(repr((node, main_child, other_child))))
        if ((main_child is not None) and (main_child['copy_number'] == node['copy_number'])):
            if (not are_all_descendants_zero(other_child)):
                raise ValueError("Invalid tree: child and complement copy_number do not add up to parent's copy_number")
            main_child['epoch_created'] = node['epoch_created']
            return main_child
        _result = node
        with open('IO_LOGGER/check_and_remove_redundant_node_io_log.txt', 'a') as f:
            f.write('Output: {0}\n'.format(repr(_result)))
        return node
    if (node is None):
        return None
    node['child'] = remove_redundant_parents(node['child'])
    node['complement'] = remove_redundant_parents(node['complement'])
    if (node['parent'] is not None):
        node = check_and_remove_redundant_node(node, node['child'], node['complement'])
        node = check_and_remove_redundant_node(node, node['complement'], node['child'])
    _result = node
    with open('IO_LOGGER/remove_redundant_parents_io_log.txt', 'a') as f:
        f.write('Output: {0}\n'.format(repr(_result)))
    return node

pass

def with_epoch_killed(tree, max_epochs):
    with open('IO_LOGGER/with_epoch_killed_io_log.txt', 'a') as f:
        f.write('Input: {0}\n'.format(repr((tree, max_epochs))))
    "\n    The function calculates the epoch_killed for each node in the given tree.\n\n    :param tree: A dictionary representing a tree node with keys:\n                 'child', 'complement', 'epoch_created', and 'epoch_killed'.\n    :return: The modified tree with 'epoch_killed' values added.\n    "
    if (tree is None):
        return None
    child = tree.get('child')
    complement = tree.get('complement')
    if ((child is not None) and (complement is not None)):
        child_epoch_created = child.get('epoch_created')
        complement_epoch_created = complement.get('epoch_created')
        if (child_epoch_created != complement_epoch_created):
            raise ValueError('Epoch created values of child and complement do not match.')
        tree['child'] = with_epoch_killed(child, max_epochs)
        tree['complement'] = with_epoch_killed(complement, max_epochs)
        tree['epoch_killed'] = child['epoch_created']
    else:
        tree['epoch_killed'] = max_epochs
    _result = tree
    with open('IO_LOGGER/with_epoch_killed_io_log.txt', 'a') as f:
        f.write('Output: {0}\n'.format(repr(_result)))
    return tree

pass

def create_truth_trees(simulated_chromosomes, max_epochs):
    with open('IO_LOGGER/create_truth_trees_io_log.txt', 'a') as f:
        f.write('Input: {0}\n'.format(repr((simulated_chromosomes, max_epochs))))
    trees = {}
    for chrom_type in simulated_chromosomes:
        sorted_list = sorted([(x['unique_identifier'], x) for x in simulated_chromosomes[chrom_type]])
        tree = {'unique_identifier': (- 1), 'parent': None, 'epoch_created': (- 1), 'paternal': None, 'child': None, 'complement': None, 'SNVs': []}
        for new_node in sorted_list:
            trees[chrom_type] = insert_node_into_truth_tree_v2(tree, new_node[1])
        for i in range(5):
            pretty_print(('##### chrom_type:' + str(chrom_type)))
        pretty_print('with nodes inserted:')
        pretty_print(trees[chrom_type])
        trees[chrom_type] = add_copynumber_to_tree(trees[chrom_type])
        pretty_print(CN_tree_from_truth_tree(trees[chrom_type]))
        pretty_print('with copynumber annotated:')
        pretty_print(trees[chrom_type])
        trees[chrom_type] = add_SNV_multiplicity_to_tree(trees[chrom_type])
        pretty_print(CN_tree_from_truth_tree(trees[chrom_type]))
        pretty_print('with SNV multiplicity:')
        pretty_print(trees[chrom_type])
        trees[chrom_type] = remove_SNVs_from_tree(trees[chrom_type])
        pretty_print(CN_tree_from_truth_tree(trees[chrom_type]))
        pretty_print('with SNVs removed from the tree:')
        pretty_print(trees[chrom_type])
        pretty_print(CN_tree_from_truth_tree(trees[chrom_type]))
        trees[chrom_type] = remove_dead_nodes(trees[chrom_type])
        pretty_print('with dead nodes removed:')
        pretty_print(trees[chrom_type])
        pretty_print(CN_tree_from_truth_tree(trees[chrom_type]))
        trees[chrom_type] = remove_redundant_parents(trees[chrom_type])
        pretty_print('with redundant nodes removed:')
        pretty_print(trees[chrom_type])
        pretty_print(CN_tree_from_truth_tree(trees[chrom_type]))
        trees[chrom_type] = with_epoch_killed(trees[chrom_type], max_epochs)
        pretty_print('with epoch killed:')
        pretty_print(trees[chrom_type])
        pretty_print(CN_tree_from_truth_tree(trees[chrom_type]))
        pretty_print(make_left_heavy(CN_tree_from_truth_tree(trees[chrom_type])))
        pretty_print(('######\n' * 5))
    _result = trees
    with open('IO_LOGGER/create_truth_trees_io_log.txt', 'a') as f:
        f.write('Output: {0}\n'.format(repr(_result)))
    return trees

def get_child_and_complement_trees(truth_tree):
    with open('IO_LOGGER/get_child_and_complement_trees_io_log.txt', 'a') as f:
        f.write('Input: {0}\n'.format(repr(truth_tree)))
    child_tree = (CN_tree_from_truth_tree(truth_tree['child']) if truth_tree['child'] else None)
    complement_tree = (CN_tree_from_truth_tree(truth_tree['complement']) if truth_tree['complement'] else None)
    _result = (child_tree, complement_tree)
    with open('IO_LOGGER/get_child_and_complement_trees_io_log.txt', 'a') as f:
        f.write('Output: {0}\n'.format(repr(_result)))
    return (child_tree, complement_tree)

def CN_tree_from_truth_tree(truth_tree):
    with open('IO_LOGGER/CN_tree_from_truth_tree_io_log.txt', 'a') as f:
        f.write('Input: {0}\n'.format(repr(truth_tree)))
    (child_tree, complement_tree) = get_child_and_complement_trees(truth_tree)
    if (child_tree and complement_tree):
        CN_tree = [truth_tree['copy_number'], child_tree, complement_tree]
    elif child_tree:
        CN_tree = [truth_tree['copy_number'], child_tree]
    elif complement_tree:
        CN_tree = [truth_tree['copy_number'], complement_tree]
    else:
        CN_tree = [truth_tree['copy_number']]
    _result = CN_tree
    with open('IO_LOGGER/CN_tree_from_truth_tree_io_log.txt', 'a') as f:
        f.write('Output: {0}\n'.format(repr(_result)))
    return CN_tree

pass

def get_left_and_right_child(tree):
    with open('IO_LOGGER/get_left_and_right_child_io_log.txt', 'a') as f:
        f.write('Input: {0}\n'.format(repr(tree)))
    left_child = tree[1]
    right_child = (tree[2] if (len(tree) == 3) else [0])
    _result = (left_child, right_child)
    with open('IO_LOGGER/get_left_and_right_child_io_log.txt', 'a') as f:
        f.write('Output: {0}\n'.format(repr(_result)))
    return (left_child, right_child)

def swap_and_make_left_heavy(tree, left_child, right_child):
    with open('IO_LOGGER/swap_and_make_left_heavy_io_log.txt', 'a') as f:
        f.write('Input: {0}\n'.format(repr((tree, left_child, right_child))))
    _result = [tree[0], make_left_heavy(right_child), make_left_heavy(left_child)]
    with open('IO_LOGGER/swap_and_make_left_heavy_io_log.txt', 'a') as f:
        f.write('Output: {0}\n'.format(repr(_result)))
    return [tree[0], make_left_heavy(right_child), make_left_heavy(left_child)]

def make_left_heavy(tree):
    with open('IO_LOGGER/make_left_heavy_io_log.txt', 'a') as f:
        f.write('Input: {0}\n'.format(repr(tree)))
    if (len(tree) == 1):
        return tree
    else:
        (left_child, right_child) = get_left_and_right_child(tree)
        if (left_child[0] < right_child[0]):
            return swap_and_make_left_heavy(tree, left_child, right_child)
        else:
            return [tree[0], make_left_heavy(left_child), make_left_heavy(right_child)]
    _result = None
    with open('IO_LOGGER/make_left_heavy_io_log.txt', 'a') as f:
        f.write('Output: {0}\n'.format(repr(_result)))

def CN_trees_from_truth_trees(truth_trees):
    with open('IO_LOGGER/CN_trees_from_truth_trees_io_log.txt', 'a') as f:
        f.write('Input: {0}\n'.format(repr(truth_trees)))
    for chrom_type in truth_trees:
        truth_trees[chrom_type] = CN_tree_from_truth_tree(truth_trees[chrom_type])
        truth_trees[chrom_type] = make_left_heavy(truth_trees[chrom_type])
    _result = truth_trees
    with open('IO_LOGGER/CN_trees_from_truth_trees_io_log.txt', 'a') as f:
        f.write('Output: {0}\n'.format(repr(_result)))
    return truth_trees

pass

pass

def count_CNs(simulated_chromosomes):
    with open('IO_LOGGER/count_CNs_io_log.txt', 'a') as f:
        f.write('Input: {0}\n'.format(repr(simulated_chromosomes)))
    observed_CNs = {}
    for chrom_type in simulated_chromosomes:
        observed_CNs[chrom_type] = [len([x for x in simulated_chromosomes[chrom_type] if ((paternal == x['paternal']) and (not x['dead']))]) for paternal in [True, False]]
    _result = observed_CNs
    with open('IO_LOGGER/count_CNs_io_log.txt', 'a') as f:
        f.write('Output: {0}\n'.format(repr(_result)))
    return observed_CNs

def count_CN_multiplicities(observed_CNs):
    with open('IO_LOGGER/count_CN_multiplicities_io_log.txt', 'a') as f:
        f.write('Input: {0}\n'.format(repr(observed_CNs)))
    multiplicities = {}
    for chrom_type in observed_CNs:
        for CN in observed_CNs[chrom_type]:
            if (CN not in multiplicities):
                multiplicities[CN] = 1
            else:
                multiplicities[CN] += 1
    _result = multiplicities
    with open('IO_LOGGER/count_CN_multiplicities_io_log.txt', 'a') as f:
        f.write('Output: {0}\n'.format(repr(_result)))
    return multiplicities

def pretty_print(message):
    with open('IO_LOGGER/pretty_print_io_log.txt', 'a') as f:
        f.write('Input: {0}\n'.format(repr(message)))
    _result = print(message)
    with open('IO_LOGGER/pretty_print_io_log.txt', 'a') as f:
        f.write('Output: {0}\n'.format(repr(_result)))
    print(message)

def random_decimal(min_val, max_val, decimal_places):
    with open('IO_LOGGER/random_decimal_io_log.txt', 'a') as f:
        f.write('Input: {0}\n'.format(repr((min_val, max_val, decimal_places))))
    _result = round(random.uniform(min_val, max_val), decimal_places)
    with open('IO_LOGGER/random_decimal_io_log.txt', 'a') as f:
        f.write('Output: {0}\n'.format(repr(_result)))
    return round(random.uniform(min_val, max_val), decimal_places)

def biased_sample(p, min_value, max_value):
    with open('IO_LOGGER/biased_sample_io_log.txt', 'a') as f:
        f.write('Input: {0}\n'.format(repr((p, min_value, max_value))))
    other_values_count = (max_value - min_value)
    if (random.random() < p):
        return min_value
    else:
        return random.randint((min_value + 1), max_value)
    _result = None
    with open('IO_LOGGER/biased_sample_io_log.txt', 'a') as f:
        f.write('Output: {0}\n'.format(repr(_result)))

def get_simulation_parameters(max_epochs):
    with open('IO_LOGGER/get_simulation_parameters_io_log.txt', 'a') as f:
        f.write('Input: {0}\n'.format(repr(max_epochs)))
    max_anue = ((max_epochs - 2) // 3)
    pre = random.randint(0, max_anue)
    mid = biased_sample(0.5, (- 1), max_anue)
    post = biased_sample(0.8, (- 1), max_anue)
    if ((mid == (- 1)) and (post >= 0)):
        (mid, post) = (post, mid)
    total_epochs = (((pre + mid) + post) + 2)
    p_up = random_decimal(0.1, 0.3, 2)
    p_down = random_decimal((p_up * 0.5), min(0.3, (p_up * 2)), 2)
    rate = 10
    _result = (pre, mid, post, p_up, p_down, rate, total_epochs)
    with open('IO_LOGGER/get_simulation_parameters_io_log.txt', 'a') as f:
        f.write('Output: {0}\n'.format(repr(_result)))
    return (pre, mid, post, p_up, p_down, rate, total_epochs)

def print_simulation_parameters(pre, mid, post, p_up, p_down, rate):
    with open('IO_LOGGER/print_simulation_parameters_io_log.txt', 'a') as f:
        f.write('Input: {0}\n'.format(repr((pre, mid, post, p_up, p_down, rate))))
    pretty_print('SIMULATION PARAMETERS ARE:')
    pretty_print(f'pre: {pre}')
    pretty_print(f'mid: {mid}')
    pretty_print(f'post: {post}')
    pretty_print(f'p_up: {p_up}')
    pretty_print(f'p_down: {p_down}')
    pretty_print(f'rate: {rate}')
    _result = pretty_print('')
    with open('IO_LOGGER/print_simulation_parameters_io_log.txt', 'a') as f:
        f.write('Output: {0}\n'.format(repr(_result)))
    pretty_print('')

def print_simulated_genome_data(simulated_chromosomes, truth_trees, CN_trees, observed_CNs, observed_CN_multiplicities):
    with open('IO_LOGGER/print_simulated_genome_data_io_log.txt', 'a') as f:
        f.write('Input: {0}\n'.format(repr((simulated_chromosomes, truth_trees, CN_trees, observed_CNs, observed_CN_multiplicities))))
    pretty_print('Simulated genome was:')
    for chrom in simulated_chromosomes:
        pretty_print(('chrom: ' + str(chrom)))
        pretty_print(simulated_chromosomes[chrom])
    pretty_print('the truth trees are:')
    for chrom_type in truth_trees:
        pretty_print(truth_trees[chrom_type])
    pretty_print('the CN simplified trees are:')
    for chrom_type in CN_trees:
        pretty_print(CN_trees[chrom_type])
    pretty_print('observed chromosomal copynumbers')
    pretty_print(observed_CNs)
    pretty_print('observed copynumber multiplicities')
    _result = pretty_print(observed_CN_multiplicities)
    with open('IO_LOGGER/print_simulated_genome_data_io_log.txt', 'a') as f:
        f.write('Output: {0}\n'.format(repr(_result)))
    pretty_print(observed_CN_multiplicities)

def save_results_to_file(test_case, simulated_chromosomes, truth_trees, CN_trees, observed_CNs, observed_CN_multiplicities, pre, mid, post, p_up, p_down, rate):
    with open('IO_LOGGER/save_results_to_file_io_log.txt', 'a') as f:
        f.write('Input: {0}\n'.format(repr((test_case, simulated_chromosomes, truth_trees, CN_trees, observed_CNs, observed_CN_multiplicities, pre, mid, post, p_up, p_down, rate))))
    file_name = f'simulation_{test_case}.txt'
    with shelve.open(file_name) as d:
        d['simulated_chromosomes'] = simulated_chromosomes
        d['truth_trees'] = truth_trees
        d['CN_trees'] = CN_trees
        d['observed_CNs'] = observed_CNs
        d['observed_CN_multiplicities'] = observed_CN_multiplicities
        d['pre'] = pre
        d['mid'] = mid
        d['post'] = post
        d['p_up'] = p_up
        d['p_down'] = p_down
        d['rate'] = rate
    _result = None
    with open('IO_LOGGER/save_results_to_file_io_log.txt', 'a') as f:
        f.write('Output: {0}\n'.format(repr(_result)))

def main():
    with open('IO_LOGGER/main_io_log.txt', 'a') as f:
        f.write('Input: {0}\n'.format(repr(())))
    pretty_print('Doing simulation')
    max_epochs = 8
    (pre, mid, post, p_up, p_down, rate, total_epochs) = get_simulation_parameters(max_epochs)
    print_simulation_parameters(pre, mid, post, p_up, p_down, rate)
    simulated_chromosomes = simulate_single_with_poisson_timestamps_names(p_up=p_up, p_down=p_down, pre=pre, mid=mid, post=post, rate=rate)
    truth_trees = create_truth_trees(simulated_chromosomes=simulated_chromosomes, max_epochs=(((pre + mid) + post) + 2))
    CN_trees = CN_trees_from_truth_trees(truth_trees)
    observed_CNs = count_CNs(simulated_chromosomes=simulated_chromosomes)
    observed_CN_multiplicities = count_CN_multiplicities(observed_CNs=observed_CNs)
    print_simulated_genome_data(simulated_chromosomes, truth_trees, CN_trees, observed_CNs, observed_CN_multiplicities)
    test_case = sys.argv[1]
    _result = save_results_to_file(test_case, simulated_chromosomes, truth_trees, CN_trees, observed_CNs, observed_CN_multiplicities, pre, mid, post, p_up, p_down, rate)
    with open('IO_LOGGER/main_io_log.txt', 'a') as f:
        f.write('Output: {0}\n'.format(repr(_result)))
    save_results_to_file(test_case, simulated_chromosomes, truth_trees, CN_trees, observed_CNs, observed_CN_multiplicities, pre, mid, post, p_up, p_down, rate)
if (__name__ == '__main__'):
    main()
