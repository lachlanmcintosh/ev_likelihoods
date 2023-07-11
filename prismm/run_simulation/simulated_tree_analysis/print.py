import logging
from typing import Dict, List, Any

def print_simulated_genome_data(
    simulated_chromosomes: Dict[str, List[Any]], 
    truth_trees: Dict[str, List[Any]], 
    CN_trees: Dict[str, List[Any]], 
    observed_CNs: Dict[str, List[Any]], 
    observed_CN_multiplicities: Dict[str, List[Any]]
) -> None:
    """
    Logs some information about simulated chromosomes and trees.

    Args:
        simulated_chromosomes: A dictionary containing lists of simulated chromosomes.
        truth_trees: A dictionary containing lists of truth trees.
        CN_trees: A dictionary containing lists of CN trees.
        observed_CNs: A dictionary containing lists of observed CNs.
        observed_CN_multiplicities: A dictionary containing lists of observed CN multiplicities.
    """
    logging.info("Simulated genome was:")
    for chrom in simulated_chromosomes:
        logging.debug("chrom: %s", str(chrom))
        logging.debug("%s", simulated_chromosomes[chrom])

    logging.debug("the truth trees are:")
    for chrom_type in truth_trees:
        logging.debug("%s", truth_trees[chrom_type])

    for chrom_type in CN_trees:
        logging.info("the CN simplified trees are:")
        logging.info("%s", CN_trees[chrom_type])

        logging.info("observed chromosomal copynumbers")
        logging.info("%s", observed_CNs[chrom_type])
