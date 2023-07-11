
def test_simulate_snvs():
    chromosomes = initialize_simulated_chromosomes()
    rate = 0.05
    epoch = 1

    initial_snv_count = 0
    for chrom_type, chrom_data in chromosomes.items():
        for entry in chrom_data:
            initial_snv_count += len(entry["SNVs"])

    final_snv_count, updated_chromosomes = simulate_snvs(chromosomes, LENGTHS, rate, epoch, initial_snv_count)
    assert final_snv_count >= initial_snv_count, "SNV count should increase after simulation"

    for chrom_type, chrom_data in updated_chromosomes.items():
        for entry in chrom_data:
            SNVS = entry["SNVs"]
            assert type(SNVS) == list
            for snv in SNVS:
                assert type(snv) == dict
                assert 'unique_identifier' in snv
                assert type(snv['unique_identifier']) == str
                assert 'epoch_created' in snv
                assert type(snv['epoch_created']) == int
                assert snv['epoch_created'] == 1


test_simulate_snvs()

from typing import Dict, List



def test_simulate_gd():
    chromosomes = initialize_simulated_chromosomes()
    epoch = 1
    initial_chrom_count = 46

    final_chrom_count, updated_chromosomes = simulate_gd(chromosomes, epoch, initial_chrom_count)
    assert final_chrom_count == initial_chrom_count * 2, "Chromosome count should double after genome doubling"

    for chrom_type, chrom_data in updated_chromosomes.items():
        for entry in chrom_data:
            if entry["epoch_created"] == epoch:
                assert entry["parent"] != -1, "Parent should not be -1 for newly created chromosomes"

    # Count the number of chromosomes in updated_chromosomes
    actual_chrom_count = sum(len(chrom_list) for chrom_list in updated_chromosomes.values())

    # Check if every chrom_type has an even number of chromosomes
    for chrom_type, chrom_list in updated_chromosomes.items():
        assert len(chrom_list) % 2 == 0, f"Chrom_type {chrom_type} should have an even number of chromosomes"

    # Check if actual chromosome count in updated_chromosomes is even
    assert actual_chrom_count % 2 == 0, "Actual chromosome count in updated_chromosomes should be an even number"

    # Check if actual chromosome count in updated_chromosomes matches the final_chrom_count
    assert actual_chrom_count == final_chrom_count, "Actual chromosome count in updated_chromosomes should match the final_chrom_count"

    # Check if half of each chromosome type have epoch = 2
    for chrom_type, chrom_list in updated_chromosomes.items():
        total_chromosomes_type = len(chrom_list)
        chromosomes_epoch_2_type = sum(chrom['epoch_created'] == 1 for chrom in chrom_list)

        assert chromosomes_epoch_2_type == total_chromosomes_type / 2, f"Half of the chromosomes of type {chrom_type} should have epoch = 2"


test_simulate_gd()

def test_simulate_anueploidy_agnostic():
    chromosomes = initialize_simulated_chromosomes()
    epoch = 1
    chrom_count = 46

    updated_chrom_count, updated_chromosomes = simulate_anueploidy_agnostic(chromosomes, epoch, chrom_count)
    assert 0 <= updated_chrom_count <= 2 * chrom_count, "Chromosome count should be within the range [0, 2 * chrom_count] after 1 round of simulation"


    # Check all newly created chromosomes have epoch == 1
    for chrom_list in updated_chromosomes.values():
        for chrom in chrom_list:
            if chrom["epoch_created"] > epoch - 1:  # for newly created chromosomes
                assert chrom["epoch_created"] == 1, "Newly created chromosomes should have epoch_created = 1"


test_simulate_anueploidy_agnostic()


def test_simulate_anueploidy():
    chromosomes = initialize_simulated_chromosomes()
    epoch = 1
    chrom_count = 46
    p_up = 0.1
    p_down = 0.1

    updated_chrom_count, updated_chromosomes = simulate_anueploidy(chromosomes, epoch, chrom_count, p_up, p_down)
    assert updated_chrom_count >= chrom_count, "Chromosome count should increase or stay the same after simulation"


    # Check all newly created chromosomes have epoch == 1
    for chrom_list in updated_chromosomes.values():
        for chrom in chrom_list:
            if chrom["epoch_created"] > epoch - 1:  # for newly created chromosomes
                assert chrom["epoch_created"] == 1, "Newly created chromosomes should have epoch_created = 1"

test_simulate_anueploidy()



def test_simulate_cancer_genome():
    p_up, p_down = 0.1, 0.1
    pre, mid, post = 2, 1, 3
    rate = 0.01
    agnostic = False

    if mid == -1:
        assert post == -1, "If mid is -1, post should also be -1"

    simulated_chromosomes = simulate_cancer_genome(
        p_up, p_down, pre, mid, post, rate, agnostic
    )

    for chrom_type in simulated_chromosomes:
        assert len(simulated_chromosomes) != 0, "There should be chromosomes in the simulation"

    if post == 0 or (mid == 0 and post == -1):
        assert simulated_chromosomes[-1] == "G", "If a genome doubling round just occurred, the last element of the simulated_chromosomes should be 'G'"

    assert check_all_chrs_are_unique(simulated_chromosomes), "All chromosomes in the simulation should be unique"
    assert check_expected_keys_in_simulated_chromosomes_present(simulated_chromosomes)


test_simulate_cancer_genome()


##### STEP 1b; from the simulated genome create a tree
#####
#####
#####
#####
#####


def test_insert_node_into_leaf():
    # Case 1: Insert a node under the child
    tree = {
        "unique_identifier": 1,
        "child": None,
        "complement": None
    }
    node = {
        "unique_identifier": 2,
        "parent": 1,
        "epoch_created": 100
    }
    expected_tree = {
        "unique_identifier": 1,
        "child": {
            "unique_identifier": 2,
            "parent": 1,
            "epoch_created": 100,
            "child": None,
            "complement": None
        },
        "complement": {
            "unique_identifier": 1,
            "child": None,
            "complement": None,
            "epoch_created": 100
        }
    }
    insert_node_into_leaf(tree, node)
    assert tree == expected_tree, f"Expected {expected_tree}, but got {tree}"


def test_insert_node_under_complement():
    # Case 1: Insert a node under the complement
    tree = {
        "unique_identifier": 1,
        "child": {
            "unique_identifier": 2,
            "parent": 1,
            "child": None,
            "complement": None
        },
        "complement": {
            "unique_identifier": 1,
            "child": None,
            "complement": None
        }
    }
    node = {
        "unique_identifier": 3,
        "parent": 1,
        "epoch_created": 1
    }
    expected_tree = {
        "unique_identifier": 1,
        "child": {
            "unique_identifier": 2,
            "parent": 1,
            "child": None,
            "complement": None
        },
        "complement": {
            "unique_identifier": 1,
            "child": {
                "unique_identifier": 3,
                "parent": 1,
                "epoch_created": 1,
                "child": None,
                "complement": None
            },
            "complement": {
                "unique_identifier": 1,
                "child": None,
                "complement": None,
                "epoch_created": 1
            }
        }
    }
    insert_node_under_complement(tree, node)
    assert tree == expected_tree, f"Expected:\n{json.dumps(expected_tree, indent=4)}\nBut got:\n{json.dumps(tree, indent=4)}"
