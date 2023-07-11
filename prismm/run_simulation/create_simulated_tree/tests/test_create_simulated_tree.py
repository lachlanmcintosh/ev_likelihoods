import pytest
from create_simulated_tree import create_truth_trees_from_simulation

def test_create_truth_trees_from_simulation():
    simulated_chromosomes = {
        'chromosome_type_1': [
            {'unique_identifier': 1, 'epoch_created': 1, 'dead': False, 'SNVs': []},
            {'unique_identifier': 2, 'epoch_created': 2, 'dead': True, 'SNVs': ['SNV_1']}
        ],
        'chromosome_type_2': [
            {'unique_identifier': 3, 'epoch_created': 1, 'dead': False, 'SNVs': ['SNV_2']},
            {'unique_identifier': 4, 'epoch_created': 2, 'dead': False, 'SNVs': []}
        ]
    }
    max_epochs = 3

    trees = create_truth_trees_from_simulation(simulated_chromosomes, max_epochs)

    assert isinstance(trees, dict)
    assert len(trees) == 2

    # Verify the truth tree for 'chromosome_type_1'
    tree_1 = trees['chromosome_type_1']
    assert tree_1 == {
        'unique_identifier': -1,
        'parent': None,
        'epoch_created': 0,
        'paternal': None,
        'child': {
            'unique_identifier': 1,
            'parent': -1,
            'epoch_created': 1,
            'paternal': True,
            'dead': False,
            'child': None,
            'complement': None,
            'SNVs': [],
            'copy_number': 1,
            'SNV_multiplicity': 0,
            'epoch_killed': 0
        },
        'complement': {
            'unique_identifier': 2,
            'parent': -1,
            'epoch_created': 2,
            'paternal': False,
            'dead': True,
            'child': None,
            'complement': None,
            'SNVs': ['SNV_1'],
            'copy_number': 1,
            'SNV_multiplicity': 0,
            'epoch_killed': 0
        },
        'copy_number': 2,
        'SNV_multiplicity': 0,
        'epoch_killed': 0
    }

    # Verify the truth tree for 'chromosome_type_2'
    tree_2 = trees['chromosome_type_2']
    assert tree_2 == {
        'unique_identifier': -1,
        'parent': None,
        'epoch_created': 0,
        'paternal': None,
        'child': {
            'unique_identifier': 3,
            'parent': -1,
            'epoch_created': 1,
            'paternal': True,
            'dead': False,
            'child': None,
            'complement': None,
            'SNVs': ['SNV_2'],
            'copy_number': 1,
            'SNV_multiplicity': 0,
            'epoch_killed': 0
        },
        'complement': {
            'unique_identifier': 4,
            'parent': -1,
            'epoch_created': 2,
            'paternal': False,
            'dead': False,
            'child': None,
            'complement': None,
            'SNVs': [],
            'copy_number': 1,
            'SNV_multiplicity': 0,
            'epoch_killed': 0
        },
        'copy_number': 2,
        'SNV_multiplicity': 0,
        'epoch_killed': 0
    }

