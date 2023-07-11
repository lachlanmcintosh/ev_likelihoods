import pytest
from mymodule import initialize_simulated_chromosomes

def test_initialize_simulated_chromosomes():
    chromosomes = initialize_simulated_chromosomes()
    assert len(chromosomes) == 23, 'Number of chromosomes should be 23'
    for chrom_type, chrom_data in chromosomes.items():
        assert len(chrom_data) == 2, f'Each chromosome type should have 2 entries. Received {len(chrom_data)} for chromosome {chrom_type}'
        for entry in chrom_data:
            expected_unique_id = chrom_type + 23 * (1 - int(entry['paternal']))
            assert entry['unique_identifier'] == expected_unique_id, f'Expected unique identifier {expected_unique_id}, but got {entry["unique_identifier"]}'
            assert entry['parent'] == -1, f'Expected parent -1, but got {entry["parent"]}'
            assert entry['epoch_created'] == 0, f'Expected epoch_created 0, but got {entry["epoch_created"]}'
            assert isinstance(entry['paternal'], bool), f'Expected paternal to be a boolean, but got {type(entry["paternal"]).__name__}'
            assert entry['SNVs'] == [], f'Expected SNVs to be an empty list, but got {entry["SNVs"]}'
            assert entry['dead'] == False, f'Expected dead to be False, but got {entry["dead"]}'
