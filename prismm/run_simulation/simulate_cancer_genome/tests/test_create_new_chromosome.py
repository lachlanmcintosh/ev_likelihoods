import pytest
import copy
from typing import Dict
from your_module_path import create_new_chromosome

def test_create_new_chromosome_typical():
    old_chromosome = {
        'unique_identifier': 0,
        'parent': -1,
        'epoch_created': 0,
        'paternal': True,
        'SNVs': [{'unique_identifier': str(i), 'epoch_created': 0} for i in range(1, 82)],
        'dead': False
    }
    chrom_count = 46
    epoch = 1
    new_chromosome = create_new_chromosome(copy.deepcopy(old_chromosome), chrom_count, epoch)
    assert new_chromosome['unique_identifier'] == 47
    assert new_chromosome['parent'] == 0
    assert new_chromosome['epoch_created'] == 1
    assert new_chromosome['paternal'] == True
    assert new_chromosome['SNVs'] == old_chromosome['SNVs']
    assert new_chromosome['dead'] == False

def test_create_new_chromosome_type_check():
    old_chromosome = {'unique_identifier': 0}
    chrom_count = "46"
    epoch = 1
    with pytest.raises(AssertionError):
        create_new_chromosome(old_chromosome, chrom_count, epoch)

def test_create_new_chromosome_empty_chromosome():
    old_chromosome = {}
    chrom_count = 46
    epoch = 1
    with pytest.raises(KeyError):
        create_new_chromosome(old_chromosome, chrom_count, epoch)
