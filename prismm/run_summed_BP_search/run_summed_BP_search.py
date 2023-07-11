import logging
from parse_arguments import parse_arguments
from IO_operations import load_results_from_file, save_likelihoods
from print_things import pretty_print_data, print_dataframes, print_path_likelihoods
from multiplicities_to_likelihoods import CN_multiplicities_to_likelihoods
from compute_likelihoods import compute_likelihoods


def generate_default_paths(max_default_path_length):
    default_paths = [str(x) for x in range(max_default_path_length)]
    default_paths += [str(x) + "G" + str(y)
                      for x in range(max_default_path_length)
                      for y in range(max_default_path_length)
                      if x + y <= max_default_path_length]

    return default_paths

def set_logging(args):
    level = logging.DEBUG if args.debug else logging.INFO
    logging.basicConfig(level=level)

def load_data_and_compute_likelihoods(args, default_paths):
    data = load_results_from_file(test_case=args.test_case, simulation_name=args.simulation_filename)
    pretty_print_data(data=data)
    likelihoods = CN_multiplicities_to_likelihoods(observed_CN_multiplicities=data['observed_CN_multiplicities'])
    computed_likelihoods = compute_likelihoods(
        likelihoods=likelihoods,
        max_number_of_solutions=args.max_number_of_solutions,
        default_paths=default_paths,
        prob_dist_filter=args.prob_dist_filter,
        path_length_diff=args.path_length_diff
    )

    return data, likelihoods, computed_likelihoods

def save_and_print_dataframes(args, data, likelihoods, computed_likelihoods, default_paths):
    marginal_likelihoods, top_likelihoods, searchable_likelihoods = computed_likelihoods
    dataframes = {
        "Marginal Likelihoods" : marginal_likelihoods,
        "Top Likelihoods" : top_likelihoods,
        "Likelihoods" : likelihoods,
        "Searchable Likelihoods" : searchable_likelihoods
    }
    print_dataframes(dataframes)
    logging.info("finished_print_dataframes")

    save_likelihoods(
        test_case=args.test_case,
        simulation_name=args.simulation_filename,
        likelihoods=likelihoods,
        marginal_likelihoods=marginal_likelihoods,
        top_likelihoods=top_likelihoods,
        searchable_likelihoods=searchable_likelihoods,
        max_number_of_solutions=args.max_number_of_solutions,
        max_default_path_length=args.max_default_path_length,
        prob_dist_filter=args.prob_dist_filter,
        path_length_diff=args.path_length_diff
    )

    print_path_likelihoods(
        likelihoods=likelihoods,
        searchable_likelihoods=searchable_likelihoods,
        marginal_likelihoods=marginal_likelihoods,
        top_likelihoods=top_likelihoods,
        default_paths=default_paths,
        data=data
    )

def main():
    args = parse_arguments()
    set_logging(args)
    default_paths = generate_default_paths(max_default_path_length=args.max_default_path_length)
    data, likelihoods, computed_likelihoods = load_data_and_compute_likelihoods(args, default_paths)
    save_and_print_dataframes(args, data, likelihoods, computed_likelihoods, default_paths)

if __name__ == "__main__":
    main()
