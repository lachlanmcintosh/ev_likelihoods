import argparse
from typing import Any
from clonal_trees.script_args import add_base_args, add_simulation_args

def parse_arguments() -> Any:
    """
    Parse command-line arguments used in the simulation process.

    :return: Namespace object with parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description="Simulate single chromosome data with Poisson timestamps."
    )

    # add the relevant arguments from clonal_trees/script_args.py
    add_base_args(parser)
    add_simulation_args(parser)

    args = parser.parse_args()

    return args


