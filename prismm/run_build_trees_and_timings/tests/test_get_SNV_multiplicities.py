
# Test case for simulated_chromosomes_to_SNV_counts
def test_simulated_chromosomes_to_SNV_counts():
    simulated_chromosomes = {
        "A": [{"dead": False, "SNVs": [{"unique_identifier": "1"}, {"unique_identifier": "2"}]},
              {"dead": True, "SNVs": [{"unique_identifier": "1"}, {"unique_identifier": "3"}]}],
        "B": [{"dead": False, "SNVs": [{"unique_identifier": "2"}, {"unique_identifier": "3"}]},
              {"dead": False, "SNVs": [{"unique_identifier": "1"}, {"unique_identifier": "4"}]}],
        "C": [{"dead": False, "SNVs": [{"unique_identifier": "2"}, {"unique_identifier": "3"}]},
              {"dead": False, "SNVs": [{"unique_identifier": "1"}, {"unique_identifier": "2"}]}]
    }

    SNV_copy_counter = simulated_chromosomes_to_SNV_counts(simulated_chromosomes)

    expected_SNV_copy_counter = {
        "A": {"1": 1, "2": 1},
        "B": {"1": 1, "2": 1, "3": 1, "4": 1},
        "C": {"1": 1, "2": 2, "3": 1}
    }

    assert SNV_copy_counter == expected_SNV_copy_counter


# Test case for SNV_counts_to_SNV_multiplicities
def test_SNV_counts_to_SNV_multiplicities():

    test_cases = [
        {
            "input": {
                "SNV_copy_counter": {
                    "A": {"1": 1, "2": 1},
                    "B": {"1": 1, "2": 1, "3": 1, "4": 1}
                }
            },
            "expected_output": {
                "A": {1: 2},
                "B": {1: 4}
            }
        },
        {
            "input": {
                "SNV_copy_counter": {
                    "A": {"1": 2},
                    "B": {"2": 2, "3": 1, "4": 1}
                }
            },
            "expected_output": {
                "A": {2: 1},
                "B": {1: 2, 2: 1}
            }
        },
        {
            "input": {
                "SNV_copy_counter": {
                    "A": {},
                    "B": {}
                }
            },
            "expected_output": {
                "A": {},
                "B": {}
            }
        },
        {
            "input": {
                "SNV_copy_counter": {
                    "A": {"1": 1},
                    "B": {"2": 2, "3": 1, "4": 1}
                }
            },
            "expected_output": {
                "A": {1: 1},
                "B": {2: 1, 1: 2}
            }
        }
    ]

    for i, test_case in enumerate(test_cases):
        input_data = test_case["input"]
        expected_output = test_case["expected_output"]
        multiplicities = SNV_counts_to_SNV_multiplicities(input_data["SNV_copy_counter"])

        assert multiplicities == expected_output, f"Test case {i} failed: Expected {expected_output}, but got {multiplicities}."



test_SNV_counts_to_SNV_multiplicities()
# Test case for SNV_counts_to_SNV_multiplicities
def test_SNV_counts_to_SNV_multiplicities():

    test_cases = [
        {
            "input": {
                "SNV_copy_counter": {
                    "A": {"1": 1, "2": 1},
                    "B": {"1": 1, "2": 1, "3": 1, "4": 1}
                }
            },
            "expected_output": {
                "A": {1: 2},
                "B": {1: 4}
            }
        },
        {
            "input": {
                "SNV_copy_counter": {
                    "A": {"1": 2},
                    "B": {"2": 2, "3": 1, "4": 1}
                }
            },
            "expected_output": {
                "A": {2: 1},
                "B": {1: 2, 2: 1}
            }
        },
        {
            "input": {
                "SNV_copy_counter": {
                    "A": {},
                    "B": {}
                }
            },
            "expected_output": {
                "A": {},
                "B": {}
            }
        },
        {
            "input": {
                "SNV_copy_counter": {
                    "A": {"1": 1},
                    "B": {"2": 2, "3": 1, "4": 1}
                }
            },
            "expected_output": {
                "A": {1: 1},
                "B": {2: 1, 1: 2}
            }
        }
    ]

    for i, test_case in enumerate(test_cases):
        input_data = test_case["input"]
        expected_output = test_case["expected_output"]
        multiplicities = SNV_counts_to_SNV_multiplicities(input_data["SNV_copy_counter"])

        assert multiplicities == expected_output, f"Test case {i} failed: Expected {expected_output}, but got {multiplicities}."



test_SNV_counts_to_SNV_multiplicities() 