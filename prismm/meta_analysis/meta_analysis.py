import pickle
import argparse
import logging
import numpy as np
import json
import glob

from clonal_trees.meta_analysis.get_best_estimates import SimulationProcessor
from clonal_trees.meta_analysis.load_simulation_results import get_filenames_sorted_by_time, process_files
from clonal_trees.meta_analysis.grid_operations import create_3x3_grid, get_correct_proportion
from clonal_trees.meta_analysis.string_distances import calculate_edit_distances, print_different_ev_strings
from clonal_trees.meta_analysis.statistics_operations import print_summary_statistics, summarize_p_similarity, summarize_plambda_similarity, aggregate_similarity_scores
from clonal_trees.meta_analysis.cutoff_methods import find_best_cutoff, find_best_accuracy_box
from clonal_trees.meta_analysis.cutoff_methods import calculate_accuracy_for_given_cutoffs


def main(simulation_filename: str, cutoffs_best_estimate: tuple, cutoffs_arbitrary_estimate: tuple) -> None:
#def main(simulation_filename: str) -> None:
    """
    Main function for processing simulation results.

    Args:
        simulation_filename (str): Name of the simulation file.
    """

    results_dict = {}

    logging.basicConfig(level=logging.INFO)
    simulation_processor = SimulationProcessor(simulation_filename)
    simulation_processor.get_best_estimates()

    sorted_filenames = get_filenames_sorted_by_time(simulation_filename)
    logging.info(sorted_filenames)

    list_of_test_cases = process_files(sorted_filenames)

    logging.info("BEST ESTIMATED")
    best_estimates = simulation_processor.get_all_best_estimates(list_of_test_cases)

    grid = create_3x3_grid(best_estimates)
    logging.info(grid)
    results_dict["grid"] = grid

    correct_proportion = get_correct_proportion(grid)
    logging.info(f"The proportion of times GD is correctly guessed is {correct_proportion}")
    results_dict["proportion_estimated"] = correct_proportion

    logging.info("Summarise levenshtein path distance:")
    lv_distances = calculate_edit_distances(best_estimates)
    logging.info("the distances are:")
    logging.info(lv_distances)
    print_summary_statistics(lv_distances, "Levenshtein Path Distance")
    results_dict["levenshtein_distances"] = lv_distances

    logging.info("Print all the pairs of ev strings that are different:")
    print_different_ev_strings(best_estimates)

    logging.info("\n\nSummarise the similarity in probability estimates")
    results_dict["p_similarity"] = summarize_p_similarity(best_estimates)

    best_estimates_zero_levenshtein = [best_estimate for best_estimate, distance in zip(best_estimates, lv_distances) if distance == 0]

    logging.info("Summarise the similarity in probability estimates for when the path difference is zero")
    results_dict["p_similarity_zero_levenshtein"] = summarize_p_similarity(best_estimates_zero_levenshtein)

    print("Summarise the similarity in the SNV rates:")
    #summarize_plambda_similarity(best_estimates)
    results_dict.update(summarize_plambda_similarity(best_estimates, "all_estimates"))

    print("Summarise the similarity in the SNV rates for when the path difference is zero:")
    results_dict.update(summarize_plambda_similarity(best_estimates_zero_levenshtein, "zero_levenshtein"))
    #summarize_plambda_similarity(best_estimates_zero_levenshtein)

    logging.info("Summarise the tree similarity statistics")
    agg_df = aggregate_similarity_scores([result['distance_dataframe'] for result in best_estimates])
    logging.info(agg_df)
    results_dict["agg_df"] = agg_df

    logging.info("Summarise the tree similarity statistics when Levenshtein distance is zero")
    agg_df_zero_levenshtein = aggregate_similarity_scores([result['distance_dataframe'] for result in best_estimates_zero_levenshtein])
    logging.info(agg_df_zero_levenshtein)
    results_dict["agg_df_zero_levenshtein"] = agg_df_zero_levenshtein

    logging.info("\nFinding the best cutoff in average Total Copy Number (TCN) for predicting 'genome_doublings':")

    cutoff1_values, cutoff2_values, accuracy_values, best_cutoffs, tcn_summaries = find_best_cutoff(best_estimates, num_cutoffs=20)
    results_dict["cutoff1_values"] = cutoff1_values
    results_dict["cutoff2_values"] = cutoff2_values
    results_dict["accuracy_values"] = accuracy_values
    results_dict["best_cutoffs"] = best_cutoffs
    results_dict["tcn_summaries"] = tcn_summaries

    np.set_printoptions(precision=3, suppress=True, threshold=np.inf, linewidth=200)
    logging.info("Cutoff 1 values: {}".format(cutoff1_values))
    logging.info("Cutoff 2 values: {}".format(cutoff2_values))
    logging.info("Accuracy values: {}".format(accuracy_values))
    logging.info("Best cutoffs: {}".format(best_cutoffs))

    cutoff1_min, cutoff1_max, cutoff2_min, cutoff2_max, best_accuracy = find_best_accuracy_box(cutoff1_values, cutoff2_values, accuracy_values)
    results_dict["best_accuracy"] = best_accuracy
    results_dict["best_accuracy_box"] = {"cutoff1_min": cutoff1_min, "cutoff1_max": cutoff1_max, "cutoff2_min": cutoff2_min, "best_accuracy": best_accuracy}
    logging.info(f"Best accuracy: {best_accuracy}")
    logging.info(f"Best accuracy box boundaries: Cutoff1_min = {cutoff1_min}, Cutoff1_max = {cutoff1_max}, Cutoff2_min = {cutoff2_min}, Cutoff2_max = {cutoff2_max}")


    logging.info("Calculating accuracy for given cutoff pairs:")
    best_estimate_accuracy = calculate_accuracy_for_given_cutoffs(best_estimates, cutoffs_best_estimate[0], cutoffs_best_estimate[1])
    arbitrary_estimate_accuracy = calculate_accuracy_for_given_cutoffs(best_estimates, cutoffs_arbitrary_estimate[0], cutoffs_arbitrary_estimate[1])

    results_dict["best_estimate_accuracy"] = best_estimate_accuracy
    results_dict["arbitrary_estimate_accuracy"] = arbitrary_estimate_accuracy

    logging.info(f"Accuracy for best cohort estimate cutoffs ({cutoffs_best_estimate[0]}, {cutoffs_best_estimate[1]}): {best_estimate_accuracy}")
    logging.info(f"Accuracy for arbitrary estimate cutoffs ({cutoffs_arbitrary_estimate[0]}, {cutoffs_arbitrary_estimate[1]}): {arbitrary_estimate_accuracy}")


    with open(f'SIMULATIONS/{simulation_filename}_summary.pkl', 'wb') as f:
        pickle.dump(results_dict, f)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--base_simulation_filename", type=str, required=True, help="Base name of the simulations")
    parser.add_argument("--cutoffs_best_estimate", type=float, nargs=2, required=True, help="Cutoff values for best cohort estimate")
    parser.add_argument("--cutoffs_arbitrary_estimate", type=float, nargs=2, required=True, help="Cutoff values for arbitrary estimate")
    args = parser.parse_args()

    main(args.base_simulation_filename, tuple(args.cutoffs_best_estimate), tuple(args.cutoffs_arbitrary_estimate))

