# File: script_args.py
import argparse
from clonal_trees.run_simulation.simulation_priors.random_number_generator import random_integer_log_scale

def add_base_args(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    parser.add_argument(
        '--test_case',
        type=int,
        default=0,
        help="The name or number of the test case you want to process."
    )
    parser.add_argument(
        '--simulation_filename',
        type=str,
        default="simulation",
        help="The name of the simulation and analysis files."
    )
    parser.add_argument(
        '--debug',
        action='store_true',
        help="Enable debug logging."
    )
    parser.add_argument(
        '--alpha',
        nargs="*",
        type=float,
        default=[20.0, 20.0, 100.0],
        help="Alpha parameters for the Dirichelet distribution for randomly sampling 'p_up' and 'p_down'. For example: --alpha 20.0 20.0 100.0. Default is [20.0, 20.0, 100.0]."
    )
    parser.add_argument(
        '--lam',
        nargs='+',
        type=float,
        default=2.0,
        help="Lam parameter for Poisson distribution for randomly sampling the mean number of events per anueploidy period (optional, default is 2)."
    )
    return parser

def add_simulation_args(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    parser.add_argument(
        "--pre",
        type=int,
        help="Pre, the parameter value for the number of anueploidy epochs per chromosome before the first round of GD if there is one."
    )
    parser.add_argument(
        "--mid",
        type=int,
        help="Mid, the parameter value for the number of anueploidy epochs per chromosome after the first round of GD and before the second if there is one."
    )
    parser.add_argument(
        "--post",
        type=int,
        help="Post."
    )
    # TODO get rid of pre, mid and post and update to just be a path argument.
    parser.add_argument(
        "--p_up",
        type=int,
        help="P_up parameter value."
    )
    parser.add_argument(
        "--p_down",
        type=int,
        help="P_down parameter value."
    )
    # TODO see that these overwrite the simulated values when specifying lam and alpha
    parser.add_argument(
        "--rate",
        type=int,
        default=random_integer_log_scale(100, 1000000),
        help="The rate at which SNVs accumulate per epoch (optional, default: random integer between 100 and 1000000 sampled on a log scale)."
    )
    parser.add_argument(
        "--max_epochs",
        type=int,
        default=8,
        help="Maximum number of anueploid/gd epochs (optional, default: 8)."
    )
    parser.add_argument(
        "--gd_probabilities",
        type=float,
        nargs='+',
        default=[0.4,0.4,0.2],
        help="Probability of gd rounds (optional, default: None)."
    )
    return parser

def bp_search_args(parser : argparse.ArgumentParser) -> argparse.ArgumentParser:
    parser.add_argument(
        '--max_number_of_solutions',
        type=int,
        default=1,
        help="The maximum number of solutions to compute (default: 20)."
    )
    parser.add_argument(
        '--max_default_path_length',
        type=int,
        default=1,
        help="The maximum length of the default paths (default: 3)."
    )
    parser.add_argument(
        '--prob_dist_filter',
        type=float,
        default=0.1,
        help="The probability distribution filter value (default: 0.1)."
    )
    parser.add_argument(
        '--path_length_diff',
        type=int,
        default=2,
        help="The path length difference value (default: 2)."
    )
    return parser

def add_build_trees_args(parser : argparse.ArgumentParser) -> argparse.ArgumentParser:
    parser.add_argument(
        '--p_window',
        type=int,
        default=0,
        help="Specify the p_window value."
    )
    parser.add_argument(
        '--plambda_window',
        type=float,
        default=0.1,
        help="Specify the plambda_window value."
    )
    parser.add_argument(
        '--tree_flexibility',
        type=int,
        default=5,
        help="Specify the tree flexibility value."
    )
    return parser

def add_analyse_trees_args(parser : argparse.ArgumentParser) -> argparse.ArgumentParser:
    return parser
