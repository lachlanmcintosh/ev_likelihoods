# File: do_all2.py
import sys
import argparse
import subprocess
import clonal_trees.script_args as script_args  
import os


def get_this_directory():
    # Get the current file path
    current_file_path = os.path.abspath(__file__)

    # Get the directory of the current file
    current_directory = os.path.dirname(current_file_path)

    return current_directory

def get_arguments_from_parser(parser):
    return [(action.option_strings, action.dest, action.default, action.type, action.choices, action.help, action.metavar)
            for action in parser._actions]

def run_script(script_name, arguments):
    print(f"DOING {script_name}")

    process = subprocess.run(["python", script_name] + arguments, env=env)
    if process.returncode != 0:
        print(f"Error in {script_name}. Exiting...")
        sys.exit(1)

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Run all processing scripts with the given test case and optional parameters."
    )
    parser.add_argument(
        "-r",
        "--run_scripts",
        nargs='+',
        type=int,
        default=[1, 2, 3, 4],
        choices=[1, 2, 3, 4],
        help="Specify which scripts to run (1, 2, 3, or 4)"
    )

    # use the functions from script_args to add the arguments
    script_args.add_base_args(parser)
    script_args.add_simulation_args(parser)
    script_args.bp_search_args(parser)
    script_args.add_build_trees_args(parser)
    script_args.add_analyse_trees_args(parser)

    return parser.parse_args()

# Import the other scripts that are required to run the pipeline
from clonal_trees.run_simulation.run_simulation import main as run_simulation
from clonal_trees.run_summed_BP_search.run_summed_BP_search import main as run_summed_BP_search
from clonal_trees.run_build_trees_and_timings.run_build_trees_and_timings import main as run_build_trees_and_timings
from clonal_trees.run_analyse_trees.run_analyse_trees import main as run_analyse_trees

def main():
    args = parse_arguments()

    debug_flag = ["--debug"] if args.debug else []
    alpha_flag = ["--alpha"] + [str(a) for a in args.alpha] if args.alpha else []
    lam_flag = ["--lam"] + [str(l) for l in args.lam] if args.lam else []

    if 1 in args.run_scripts:
        run_simulation.main([
            "--test_case", str(args.test_case),
            "--simulation_filename", args.simulation_filename,
            "--rate", str(10000), # TODO remove this once debugging is done
            *debug_flag, *alpha_flag, *lam_flag
        ])

    if 2 in args.run_scripts:
        run_summed_BP_search.main([
            "--test_case", str(args.test_case),
            "--simulation_filename", args.simulation_filename,
            *debug_flag,
            "--max_number_of_solutions", str(args.max_number_of_solutions),
            "--max_default_path_length", str(args.max_default_path_length),
            "--prob_dist_filter", str(args.prob_dist_filter),
            "--path_length_diff", str(args.path_length_diff)
        ])

    if 3 in args.run_scripts:
        run_build_trees_and_timings.main([
            "--test_case", str(args.test_case),
            "--simulation_filename", args.simulation_filename,
            "--p_window", str(args.p_window),
            "--plambda_window", str(args.plambda_window),
            "--tree_flexibility", str(args.tree_flexibility),
            *debug_flag,
            *alpha_flag
        ])

    if 4 in args.run_scripts:
        run_analyse_trees.main([
            "--test_case", str(args.test_case),
            "--simulation_filename", args.simulation_filename,
            *debug_flag
        ])

if __name__ == "__main__":
    main()
