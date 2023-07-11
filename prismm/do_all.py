import sys
import argparse
import subprocess


def run_script(script_name, arguments):
    print(f"DOING {script_name}")
    process = subprocess.run(["python", script_name] + arguments)
    if process.returncode != 0:
        print(f"Error in {script_name}. Exiting...")
        sys.exit(1)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Run all processing scripts with the given test case and optional parameters."
    )
    parser.add_argument(
        "-t", 
        "--test_case", 
        type=int, 
        default=0, 
        help="The name or number of the test case you want to process."
    )
    parser.add_argument(
        "-f", 
        "--simulation_filename", 
        type=str, 
        default="simulation", 
        help="The name of the simulation and analysis files."
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
    parser.add_argument(
        "-b", 
        "--debug", 
        action='store_true', 
        help="Enable debug logging."
    )
    parser.add_argument(
        "-s", 
        "--max_number_of_solutions", 
        type=int, 
        default=1, 
        help="The maximum number of solutions to compute (default: 20)."
    )
    parser.add_argument(
        "-d", 
        "--max_default_path_length", 
        type=int, 
        default=1, 
        help="The maximum length of the default paths (default: 3)."
    )
    parser.add_argument(
        "-p", 
        "--prob_dist_filter", 
        type=float, 
        default=0.1, 
        help="The probability distribution filter value (default: 0.1)."
    )
    parser.add_argument(
        "-l", 
        "--path_length_diff", 
        type=int, 
        default=2, 
        help="The path length difference value (default: 2)."
    )
    parser.add_argument(
        "-w", 
        "--p_window", 
        type=int, 
        default=0, 
        help="Specify the p_window value."
    )
    parser.add_argument(
        "-o", 
        "--plambda_window", 
        type=float, 
        default=0.1, 
        help="Specify the plambda_window value."
    )
    parser.add_argument(
        "-x", 
        "--tree_flexibility", 
        type=int, 
        default=5, 
        help="Specify the tree flexibility value."
    )
    parser.add_argument(
        "--alpha", 
        nargs="*", 
        type=float, 
        default=[20.0, 20.0, 100.0], 
        help="Alpha values for the simulation as a list of floats. For example: --alpha 20.0 20.0 100.0. Default is [20.0, 20.0, 100.0]."
    )
    parser.add_argument(
        "-m", 
        "--lam", 
        nargs='+', 
        type=float, 
        help="Lam parameter for Dirichlet distribution (optional, default is [2, 2, 10])."
    )
    return parser.parse_args()


def main():
    args = parse_arguments()

    debug_flag = ["--debug"] + [args.debug] if args.debug else []  
    alpha_flag = ["--alpha"] + [str(a) for a in args.alpha] if args.alpha else []
    lam_flag = ["--lam"] + [str(l) for l in args.lam] if args.lam else []

    if 1 in args.run_scripts:
        run_script("run_simulation_old.py", [
            "--test_case", str(args.test_case),
            "--simulation_filename", args.simulation_filename,
            *debug_flag, *alpha_flag, *lam_flag
        ])

    if 2 in args.run_scripts:
        run_script("run_summed_BP_search.py", [
            "--test_case", str(args.test_case),
            "--simulation_filename", args.simulation_filename,
            *debug_flag,
            "--max_number_of_solutions", str(args.max_number_of_solutions),
            "--max_default_path_length", str(args.max_default_path_length),
            "--prob_dist_filter", str(args.prob_dist_filter),
            "--path_length_diff", str(args.path_length_diff)
        ])

    if 3 in args.run_scripts:
        run_script("run_build_trees_and_timings.py", [
            "--test_case", str(args.test_case),
            "--simulation_filename", args.simulation_filename,
            "--p_window", str(args.p_window),
            "--plambda_window", str(args.plambda_window),
            "--tree_flexibility", str(args.tree_flexibility),
            *debug_flag,
            *alpha_flag
        ])

    if 4 in args.run_scripts:
        run_script("run_analyse_trees.py", [
            "--test_case", str(args.test_case),
            "--simulation_filename", args.simulation_filename,
            *debug_flag
        ])

    #if 5 in args.run_scripts:
    #    run_script("do_meta_analysis.py", [
    #        "--simulation_filename", args.simulation_filename,
    #        *debug_flag
    #    ])


if __name__ == "__main__":
    main()

