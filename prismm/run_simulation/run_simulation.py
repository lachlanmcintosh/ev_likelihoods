import pickle
import copy
import logging
from typing import Dict, List, Any, Tuple
from pathlib import Path

from clonal_trees.run_simulation.simulate_cancer_genome.simulate_cancer_genome import simulate_cancer_genome
from clonal_trees.run_simulation.simulated_tree_analysis.count_multiplicities import count_copy_numbers, count_copy_number_multiplicities
from clonal_trees.run_simulation.create_simulated_tree.create_simulated_tree import create_truth_trees_from_simulation
from clonal_trees.run_simulation.create_simulated_tree.dict_tree_to_CN_tree import convert_truth_trees_to_CN_trees
from clonal_trees.run_simulation.simulated_tree_analysis.print import print_simulated_genome_data
from clonal_trees.run_simulation.simulation_priors.simulate_prior import simulate_parameters_not_given_as_arguments
from clonal_trees.run_simulation.parse_arguments import parse_arguments

# Define the simulation file folder as a Path object for easier file handling
SIMULATIONS_FILE_FOLDER = Path("clonal_trees/SIMULATIONS/")

def simulate_and_analyze_genome(args) -> Tuple[Dict, List, List, Dict, Dict]:
    """
    Simulate the genome and perform basic analysis.

    This function simulates a cancer genome based on the provided arguments, creates truth trees from the simulation,
    converts these truth trees to CN trees, and counts the copy numbers and their multiplicities.

    Args:
        args: An Arguments object containing the simulation parameters.

    Returns:
        A tuple containing the simulated chromosomes, truth trees, CN trees, observed CNs, and observed CN multiplicities.
    """
    # Simulate the genome
    print(args)
    simulated_chromosomes = simulate_cancer_genome(
        p_up=args.p_up,
        p_down=args.p_down,
        pre=args.pre,
        mid=args.mid,
        post=args.post,
        rate=args.rate
    )

    # Perform basic analysis of the simulated genome
    truth_trees = create_truth_trees_from_simulation(
        simulated_chromosomes=simulated_chromosomes,
        max_epochs=args.pre + args.mid + args.post + 2
    )

    CN_trees = convert_truth_trees_to_CN_trees(truth_trees=copy.deepcopy(truth_trees))
    observed_CNs = count_copy_numbers(simulated_chromosomes=simulated_chromosomes)
    observed_CN_multiplicities = count_copy_number_multiplicities(observed_copy_numbers=observed_CNs)

    return simulated_chromosomes, truth_trees, CN_trees, observed_CNs, observed_CN_multiplicities

# Define the simulation file folder as a Path object for easier file handling
SIMULATIONS_FILE_FOLDER = Path("clonal_trees/SIMULATIONS/")


def save_results_to_file(test_case: str, simulation_filename: str, simulated_chromosomes: Dict,
                         truth_trees: List, CN_trees: List, observed_CNs: Dict,
                         observed_CN_multiplicities: Dict, pre: int, mid: int,
                         post: int, p_up: float, p_down: float, rate: float) -> None:
    """
    Function to save the simulation results to a pickle file.
    
    Arguments:
    test_case -- string representing the test case
    simulation_filename -- string representing the filename for the simulation
    simulated_chromosomes -- dictionary representing the simulated chromosomes
    truth_trees -- list representing the truth trees
    CN_trees -- list representing the CN trees
    observed_CNs -- dictionary representing the observed CNs
    observed_CN_multiplicities -- dictionary representing the observed CN multiplicities
    pre -- integer representing the pre epoch
    mid -- integer representing the mid epoch
    post -- integer representing the post epoch
    p_up -- float representing the p_up
    p_down -- float representing the p_down
    rate -- float representing the rate
    """
    
    file_name = SIMULATIONS_FILE_FOLDER / f'{simulation_filename}_{test_case}.pickle'
    with file_name.open('wb') as f:
        pickle.dump({
            'simulated_chromosomes': simulated_chromosomes,
            'truth_trees': truth_trees,
            'CN_trees': CN_trees,
            'observed_CNs': observed_CNs,
            'observed_CN_multiplicities': observed_CN_multiplicities,
            'pre': pre,
            'mid': mid,
            'post': post,
            'p_up': p_up,
            'p_down': p_down,
            'rate': rate
        }, f)
    logging.info(f"Results saved to file {file_name}")

def main() -> None:
    """
    Main function to run the simulation, analyze the results, and save them.

    The function performs the following steps:
    1. Parse command-line arguments.
    2. Simulate parameters that are not provided as arguments.
    3. Simulate a genome and perform basic analysis.
    4. Print the simulated genome data.
    5. Save the results to a file.
    """
    
    # Get simulation parameters
    args = parse_arguments()

    args = simulate_parameters_not_given_as_arguments(args)

    simulated_chromosomes, truth_trees, CN_trees, observed_CNs, observed_CN_multiplicities = simulate_and_analyze_genome(args)

    print_simulated_genome_data(
        simulated_chromosomes=simulated_chromosomes,
        truth_trees=truth_trees,
        CN_trees=CN_trees,
        observed_CNs=observed_CNs,
        observed_CN_multiplicities=observed_CN_multiplicities
    )

    # Save the results
    save_results_to_file(
        test_case=args.test_case,
        simulation_filename=args.simulation_filename,
        simulated_chromosomes=simulated_chromosomes,
        truth_trees=truth_trees,
        CN_trees=CN_trees,
        observed_CNs=observed_CNs,
        observed_CN_multiplicities=observed_CN_multiplicities,
        pre=args.pre,
        mid=args.mid,
        post=args.post,
        p_up=args.p_up,
        p_down=args.p_down,
        rate=args.rate
    )


if __name__ == "__main__":
    main()
