from typing import List

from clonal_trees.run_simulation.simulate_cancer_genome.initialize_simulated_chromosomes import initialize_simulated_chromosomes
from clonal_trees.run_simulation.simulate_cancer_genome.simulate_SNVs import simulate_snvs
from clonal_trees.run_simulation.simulate_cancer_genome.simulation_checks import check_all_chrs_are_unique, check_simulated_chromosomes
from clonal_trees.run_simulation.simulate_cancer_genome.simulate_GD import simulate_gd
from clonal_trees.run_simulation.simulate_cancer_genome.simulate_anueploidy_agnostic import simulate_anueploidy_agnostic
from clonal_trees.run_simulation.simulate_cancer_genome.simulate_anueploidy_model import simulate_anueploidy
from clonal_trees.utils.LENGTHS import LENGTHS
from clonal_trees.utils.get_ev_string import get_ev_string

def simulate_cancer_genome(p_up: float, p_down: float, pre: int, mid: int, post: int, rate: float, agnostic: bool=False) -> List:
    """
    Simulate a genome based on given parameters.

    :param p_up: Probability of upward mutation.
    :param p_down: Probability of downward mutation.
    :param pre: Pre-mutation period.
    :param mid: Mid-mutation period.
    :param post: Post-mutation period.
    :param rate: Mutation rate.
    :param agnostic: Indicates if the simulation is agnostic or not.
    :return: Simulated chromosomes.
    """
    # Count of unique SNVs
    snv_count = 0

    # Initialize chromosomes
    simulated_chromosomes = initialize_simulated_chromosomes()

    # Count of unique chromosomes
    chrom_count = 46  # in a human cell

    # Total epochs count
    ev_sequence = get_ev_string(pre, mid, post)

    if len(ev_sequence) == 0:
        return simulated_chromosomes

    for epoch, epoch_type in enumerate(ev_sequence):
        # Simulate SNVs
        # note that SNVs are created in this epoch, whilst we "time" the copy number changes to be in the next epoch when they are created.
        snv_count, simulated_chromosomes = simulate_snvs(simulated_chromosomes, LENGTHS, rate, epoch, snv_count)
        check_all_chrs_are_unique(simulated_chromosomes)

        if (mid != -1 and epoch == pre) or (post != -1 and epoch == pre + 1 + mid):
            assert epoch_type == "G"
            chrom_count, simulated_chromosomes = simulate_gd(simulated_chromosomes, epoch+1, chrom_count)
            check_all_chrs_are_unique(simulated_chromosomes)
        else:
            assert epoch_type == "A"
            if agnostic:
                chrom_count, simulated_chromosomes = simulate_anueploidy_agnostic(simulated_chromosomes, epoch+1, chrom_count)
                check_all_chrs_are_unique(simulated_chromosomes)
            else:
                chrom_count, simulated_chromosomes = simulate_anueploidy(simulated_chromosomes, epoch+1, chrom_count, p_up, p_down)
                check_all_chrs_are_unique(simulated_chromosomes)

    assert pre + mid + post + 2 == len(ev_sequence)
    check_simulated_chromosomes(simulated_chromosomes, pre, mid, post, ev_sequence)

    return simulated_chromosomes
