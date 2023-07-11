

def test_count_paternity():
    Input1 = ([
        {
            'unique_identifier': 22,
            'parent': -1,
            'epoch_created': 0,
            'paternal': True,
            'SNVs': [],
            'dead': False
        },
        {
            'unique_identifier': 45,
            'parent': -1,
            'epoch_created': 0,
            'paternal': False,
            'SNVs': [],
            'dead': False
        },
        {
            'unique_identifier': 91,
            'parent': 22,
            'epoch_created': 0,
            'paternal': True,
            'SNVs': [],
            'dead': False
        },
        {
            'unique_identifier': 92,
            'parent': 45,
            'epoch_created': 0,
            'paternal': False,
            'SNVs': [],
            'dead': False
        }
    ], False)
    Output1 = 2

    Input2 = ([
        {
            'unique_identifier': 19,
            'parent': -1,
            'epoch_created': 0,
            'paternal': True,
            'SNVs': [{'unique_identifier': '31', 'epoch_created': 2}],
            'dead': False
        },
        {
            'unique_identifier': 42,
            'parent': -1,
            'epoch_created': 0,
            'paternal': False,
            'SNVs': [],
            'dead': True
        },
        {
            'unique_identifier': 85,
            'parent': 19,
            'epoch_created': 0,
            'paternal': True,
            'SNVs': [],
            'dead': False
        },
        {
            'unique_identifier': 106,
            'parent': 86,
            'epoch_created': 1,
            'paternal': False,
            'SNVs': [],
            'dead': False
        },
        {
            'unique_identifier': 86,
            'parent': 42,
            'epoch_created': 0,
            'paternal': False,
            'SNVs': [],
            'dead': False
        },
        {
            'unique_identifier': 189,
            'parent': 19,
            'epoch_created': 2,
            'paternal': True,
            'SNVs': [{'unique_identifier': '31', 'epoch_created': 2}],
            'dead': False
        },
        {
            'unique_identifier': 190,
            'parent': 85,
            'epoch_created': 2,
            'paternal': True,
            'SNVs': [],
            'dead': False
        },
        {
            'unique_identifier': 191,
            'parent': 106,
            'epoch_created': 2,
            'paternal': False,
            'SNVs': [],
            'dead': False
        },
        {
            'unique_identifier': 192,
            'parent': 86,
            'epoch_created': 2,
            'paternal': False,
            'SNVs': [],
            'dead': False
        }
    ], True)
    Output2 = 4

    assert count_paternity(Input1[0], Input1[1]) == Output1
    assert count_paternity(Input2[0], Input2[1]) == Output2


def test_check_all_chrs_are_unique():
    """
    Test the check_all_chrs_are_unique function.
    """
    unique_chrs = {
        "type_A": [{"unique_identifier": 1}, {"unique_identifier": 2}],
        "type_B": [{"unique_identifier": 3}, {"unique_identifier": 4}],
    }
    assert check_all_chrs_are_unique(unique_chrs)

    non_unique_chrs = {
        "type_A": [{"unique_identifier": 1}, {"unique_identifier": 2}],
        "type_B": [{"unique_identifier": 2}, {"unique_identifier": 4}],
    }
    with pytest.raises(ValueError):
        check_all_chrs_are_unique(non_unique_chrs)

    # More complex test cases
    non_unique_chrs_complex = {
        "type_A": [{"unique_identifier": i} for i in range(10)],
        "type_B": [{"unique_identifier": i} for i in range(5, 15)],
    }
    with pytest.raises(ValueError):
        check_all_chrs_are_unique(non_unique_chrs_complex)

    non_unique_chrs_single_type = {
        "type_A": [{"unique_identifier": 1}, {"unique_identifier": 1}],
    }
    with pytest.raises(ValueError):
        check_all_chrs_are_unique(non_unique_chrs_single_type)


test_check_all_chrs_are_unique()


def test_check_expected_keys_in_simulated_chromosomes_present():
    valid_chrs = {
        "type_A": [
            {
                "SNVs": [],
                "paternal": "A",
                "epoch_created": 0,
                "parent": None,
                "unique_identifier": 1,
                "dead": False,
            }
        ],
    }
    assert check_expected_keys_in_simulated_chromosomes_present(valid_chrs)

    invalid_chrs = {
        "type_A": [
            {
                "SNVs": [],
                "paternal": "A",
                "epoch_created": 0,
                "unique_identifier": 1,
            }
        ],
    }
    assert not check_expected_keys_in_simulated_chromosomes_present(invalid_chrs)

test_check_expected_keys_in_simulated_chromosomes_present()


def test_initialize_simulated_chromosomes():
    chromosomes = initialize_simulated_chromosomes()

    assert len(chromosomes) == 23, "Number of chromosomes should be 23"

    for chrom_type, chrom_data in chromosomes.items():
        assert len(chrom_data) == 2, "Each chromosome type should have 2 entries"

        for entry in chrom_data:
            assert entry["unique_identifier"] == chrom_type + 23*(1-int(entry["paternal"])), f"Unique identifier should be {chrom_type} + 0 for paternal, and {chrom_type} + 23 for maternal (entry: {entry})"
            assert entry["parent"] == -1, f"Parent of root node should be -1 (entry: {entry})"
            assert entry["epoch_created"] == 0, f"Epoch created should be 0 (entry: {entry})"
            assert isinstance(entry["paternal"], bool), f"Paternal attribute should be a boolean value (entry: {entry})"
            assert entry["SNVs"] == [], f"SNVs should be initialized as an empty list (entry: {entry})"
            assert entry["dead"] == False, f"Dead attribute should be False initially (entry: {entry})"


test_initialize_simulated_chromosomes()



