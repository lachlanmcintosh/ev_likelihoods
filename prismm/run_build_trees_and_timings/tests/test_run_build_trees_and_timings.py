def test_sum_SNV_counts():
    # Test case 1
    observed_SNV_multiplicities1 = {
        "chr1": {0: 5, 1: 10},
        "chr2": {1: 7, 2: 3},
    }
    expected_total1 = 25
    assert sum_SNV_counts(observed_SNV_multiplicities1) == expected_total1

    # Test case 2
    observed_SNV_multiplicities2 = {
        "chr1": {0: 3, 1: 4},
        "chr2": {2: 1, 3: 2},
    }
    expected_total2 = 10
    assert sum_SNV_counts(observed_SNV_multiplicities2) == expected_total2



test_sum_SNV_counts()


def test_sum_chrom_multiplicities():
    # Test case 1
    observed_CN_multiplicities1 = {"chr1": 5, "chr2": 10}
    expected_total1 = 15
    assert sum_chrom_multiplicities(observed_CN_multiplicities1) == expected_total1

    # Test case 2
    observed_CN_multiplicities2 = {"chr1": 3, "chr2": 7}
    expected_total2 = 10
    assert sum_chrom_multiplicities(observed_CN_multiplicities2) == expected_total2



test_sum_chrom_multiplicities()


def test_sum_observed_CNs():
    test_cases = [
        {
            "input": {0: [1, 0], 1: [2, 0], 2: [2, 2]},
            "expected_output": 7,
        },
        {
            "input": {3: [3, 0], 4: [1, 0], 5: [1, 0]},
            "expected_output": 5,
        },
        {
            "input": {6: [1, 0], 7: [5, 1], 8: [3, 3]},
            "expected_output": 13,
        },
        {
            "input": {9: [2, 1], 10: [1, 0], 11: [1, 0]},
            "expected_output": 5,
        },
    ]

    for test_case in test_cases:
        input_dictionary = test_case["input"]
        expected_output = test_case["expected_output"]
        result = sum_observed_CNs(input_dictionary)
        assert result == expected_output, f"For {input_dictionary}, expected {expected_output} but got {result}"


test_sum_observed_CNs()
