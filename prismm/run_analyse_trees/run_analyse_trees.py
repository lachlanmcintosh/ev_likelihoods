import argparse
from typing import List


def main(test_case: str, simulation_filename: str) -> None:
    """
    The main function that processes a simulation test case and analyzes estimated trees.

    Args:
        test_case (str): A string specifying which test case to process.
        simulation_filename (str): The filename containing the simulation data.

    """
    SS = load_simulation_data(test_case, simulation_filename)
    SS["solutions"] = sort_simulation_results_by_likelihood(SS["solutions"])
    get_all_trees_ready_for_comparison(SS)
    print_solution_results(SS)

    # List of node distance functions
    node_distance_functions = [ZeroDistance, AbsoluteDifference, SquaredDistance, InfiniteDistance]

    # Call the function
    compute_solution_similarity(SS, node_distance_functions)
    add_relative_timing_comparison(SS)
    print_similarity_results(SS)

    process_further(SS)

    save_simulation_data(test_case, simulation_filename, SS)

    
    print_dicts_memory(SS)

if __name__ == "__main__":
    args = parse_arguments()
    test_case = args.test_case
    simulation_filename = args.simulation_filename
    main(test_case, simulation_filename)


