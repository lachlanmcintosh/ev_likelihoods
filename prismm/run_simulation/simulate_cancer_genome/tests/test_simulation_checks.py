import pytest
from typing import List, Dict, Any

from simulation_checks import count_paternity

chromosome_data = [
    (
        [
            {'unique_identifier': 22, 'parent': -1, 'epoch_created': 0, 'paternal': True, 'SNVs': [], 'dead': False},
            {'unique_identifier': 45, 'parent': -1, 'epoch_created': 0, 'paternal': False, 'SNVs': [], 'dead': False},
            {'unique_identifier': 91, 'parent': 22, 'epoch_created': 0, 'paternal': True, 'SNVs': [], 'dead': False},
            {'unique_identifier': 92, 'parent': 45, 'epoch_created': 0, 'paternal': False, 'SNVs': [], 'dead': False}
        ],
        True,
        2
    ),
    (
        [
            {'unique_identifier': 19, 'parent': -1, 'epoch_created': 0, 'paternal': True, 'SNVs': [{'unique_identifier': '31', 'epoch_created': 2}], 'dead': False},
            {'unique_identifier': 42, 'parent': -1, 'epoch_created': 0, 'paternal': False, 'SNVs': [], 'dead': True},
            {'unique_identifier': 85, 'parent': 19, 'epoch_created': 0, 'paternal': True, 'SNVs': [], 'dead': False},
            {'unique_identifier': 106, 'parent': 86, 'epoch_created': 1, 'paternal': False, 'SNVs': [], 'dead': False},
            {'unique_identifier': 86, 'parent': 42, 'epoch_created': 0, 'paternal': False, 'SNVs': [], 'dead': False},
            {'unique_identifier': 189, 'parent': 19, 'epoch_created': 2, 'paternal': True, 'SNVs': [{'unique_identifier': '31', 'epoch_created': 2}], 'dead': False},
            {'unique_identifier': 190, 'parent': 85, 'epoch_created': 2, 'paternal': True, 'SNVs': [], 'dead': False},
            {'unique_identifier': 191, 'parent': 106, 'epoch_created': 2, 'paternal': False, 'SNVs': [], 'dead': False},
            {'unique_identifier': 192, 'parent': 86, 'epoch_created': 2, 'paternal': False, 'SNVs': [], 'dead': False}
        ],
        False,
        4
    ),
]

@pytest.mark.parametrize('chromosomes, paternal, expected', chromosome_data)
def test_count_paternity_valid(chromosomes, paternal, expected):
    assert count_paternity(chromosomes, paternal) == expected

def test_count_paternity_missing_paternal():
    with pytest.raises(ValueError, match="Chromosome is missing 'paternal' key"):
        count_paternity([{'unique_identifier': 22, 'parent': -1, 'epoch_created': 0, 'SNVs': [], 'dead': False}], True)

def test_count_paternity_missing_dead():
    with pytest.raises(ValueError, match="Chromosome is missing 'dead' key"):
        count_paternity([{'unique_identifier': 22, 'parent': -1, 'epoch_created': 0, 'paternal': True, 'SNVs': []}], True)

from mymodule import check_all_chrs_are_unique

def test_check_all_chrs_are_unique():
    unique_chrs = {
        'type_A': [{'unique_identifier': 1}, {'unique_identifier': 2}],
        'type_B': [{'unique_identifier': 3}, {'unique_identifier': 4}]
    }
    assert check_all_chrs_are_unique(unique_chrs)

def test_check_all_chrs_are_unique_empty_dictionary():
    assert check_all_chrs_are_unique({})

@pytest.mark.parametrize(
    'non_unique_chrs, exception',
    [
        (
            {
                'type_A': [{'unique_identifier': 1}, {'unique_identifier': 2}],
                'type_B': [{'unique_identifier': 2}, {'unique_identifier': 4}]
            },
            ValueError
        ),
        (
            {
                'type_A': [{'unique_identifier': i} for i in range(10)],
                'type_B': [{'unique_identifier': i} for i in range(5, 15)]
            },
            ValueError
        ),
        (
            {
                'type_A': [{'unique_identifier': 1}, {'unique_identifier': 1}]
            },
            ValueError
        )
    ]
)
def test_check_all_chrs_are_unique_with_errors(non_unique_chrs, exception):
    with pytest.raises(exception):
        check_all_chrs_are_unique(non_unique_chrs)

def test_check_all_chrs_are_unique_missing_key():
    with pytest.raises(KeyError):
        check_all_chrs_are_unique({'type_A': [{'other_key': 1}]})

from your_module import check_expected_keys_in_simulated_chromosomes_present

@pytest.mark.parametrize(
    "simulated_chromosomes, expected",
    [
        (
            {
                'type_A': [
                    {'SNVs': [], 'paternal': 'A', 'epoch_created': 0, 'parent': None, 'unique_identifier': 1, 'dead': False}
                ]
            },
            True
        ),
        (
            {
                'type_A': [
                    {'SNVs': [], 'paternal': 'A', 'epoch_created': 0, 'unique_identifier': 1}
                ]
            },
            False
        ),
        (
            {
                'type_A': [
                    {'SNVs': [], 'paternal': 'A', 'epoch_created': 0, 'parent': None, 'unique_identifier': 1, 'dead': False, 'extra_key': 'unexpected'}
                ]
            },
            True
        ),
        (
            {}, 
            True
        ),
        (
            {
                'type_A': []
            }, 
            True
        ),
    ],
)
def test_check_expected_keys_in_simulated_chromosomes_present(simulated_chromosomes: Dict[str, List[Dict[str, str]]], expected: bool):
    assert check_expected_keys_in_simulated_chromosomes_present(simulated_chromosomes) == expected

from your_module import check_simulated_chromosomes

def test_check_simulated_chromosomes():
    simulated_chromosomes = {
        'type_A': [{'SNVs': [], 'paternal': 'A', 'epoch_created': 0, 'parent': None, 'unique_identifier': 1, 'dead': False}]
    }
    pre = 1
    mid = 1
    post = 1
    ev_sequence = 'ABCDEFG'
    check_simulated_chromosomes(simulated_chromosomes, pre, mid, post, ev_sequence)

@pytest.mark.parametrize(
    'simulated_chromosomes, pre, mid, post, ev_sequence, exception',
    [
        (
            {},
            3,
            1,
            1,
            'ABCDEFG',
            AssertionError
        ),
        (
            {
                'type_A': [{'SNVs': [], 'paternal': 'A', 'epoch_created': 0, 'parent': None, 'unique_identifier': 1, 'dead': False},
                {'SNVs': [], 'paternal': 'A', 'epoch_created': 0, 'parent': None, 'unique_identifier': 1, 'dead': False}]
            },
            3,
            1,
            1,
            'ABCDEFG',
            AssertionError
        ),
        (
            {
                'type_A': [{'SNVs': [], 'paternal': 'A', 'epoch_created': 0, 'unique_identifier': 1}]
            },
            3,
            1,
            1,
            'ABCDEFG',
            AssertionError
        ),
        (
            {
                'type_A': [{'SNVs': [], 'paternal': 'A', 'epoch_created': 0, 'parent': None, 'unique_identifier': 1, 'dead': False}]
            },
            3,
            1,
            1,
            'ABCD',
            AssertionError
        )
    ]
)
def test_check_simulated_chromosomes_with_errors(simulated_chromosomes, pre, mid, post, ev_sequence, exception):
    with pytest.raises(exception):
        check_simulated_chromosomes(simulated_chromosomes, pre, mid, post, ev_sequence)
