import argparse
import clonal_trees.script_args as script_args

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Process and compute likelihoods for a given test case."
    )

    # add the relevant arguments from clonal_trees/script_args.py
    script_args.add_base_args(parser)
    script_args.add_BP_search_args(parser)

    return parser.parse_args()
