import logging
from typing import Union

def print_simulation_parameters(
    pre: Union[int, float], 
    mid: Union[int, float], 
    post: Union[int, float], 
    p_up: Union[int, float], 
    p_down: Union[int, float], 
    rate: Union[int, float]
) -> None:
    """
    Logs some information about simulation parameters.

    Args:
        pre: The 'pre' parameter for the simulation.
        mid: The 'mid' parameter for the simulation.
        post: The 'post' parameter for the simulation.
        p_up: The 'p_up' parameter for the simulation.
        p_down: The 'p_down' parameter for the simulation.
        rate: The 'rate' parameter for the simulation.
    """
    logging.info("SIMULATION PARAMETERS ARE:")
    logging.info("pre: %s", pre)
    logging.info("mid: %s", mid)
    logging.info("post: %s", post)
    logging.info("p_up: %s", p_up)
    logging.info("p_down: %s", p_down)
    logging.info("rate: %s", rate)
    logging.info("")
