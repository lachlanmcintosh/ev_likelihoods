import argparse
from script_args import add_simulation_args, add_analyse_trees_args

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Run all processing scripts with the given test case and optional parameters."
    )

    # add the relevant arguments from clonal_trees/script_args.py
    script_args.add_base_args(parser)
    script_args.add_analyse_trees_args(parser)    

    return parser.parse_args()