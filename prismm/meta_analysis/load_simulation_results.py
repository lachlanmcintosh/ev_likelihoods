import os
import glob
import pickle
import re
import copy
import logging
from typing import List, Optional, Dict, Any, Tuple

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class SimulationResult:
    def __init__(self, filename: str):
        self.filename = filename
        self.data = self._load_file()

    def _load_file(self) -> Optional[Dict[str, Any]]:
        """
        Load the simulation result file.

        Args:
            filename: The name of the file.

        Returns:
            The data in the file as a Python object, or None if the file couldn't be loaded.
        """
        try:
            with open(self.filename, 'rb') as f:
                return pickle.load(f)
        except EOFError:
            logger.error(f"EOFError: {self.filename} is incomplete. The file might still be generating.")
        except Exception as e:
            logger.error(f"An error occurred while loading {self.filename}: {e}")
        return None

    def _get_test_case(self) -> int:
        """
        Extract the test case number from the filename.

        Returns:
            The test case number.
        """
        return int(re.search(r'(\d+)', self.filename).group())

    def _add_test_case_to_results(self) -> None:
        """
        Add the test case number to all the results in the simulation.
        """
        if self.data and "solutions" in self.data:
            for result in self.data["solutions"]:
                result['test_case'] = self._get_test_case()

    def _print_result_info(self, result: Dict[str, Any]) -> None:
        """
        Print the truth and estimated values for each field in a result.

        Args:
            result: A dictionary representing the result of a simulation.
        """
        fields = ["pre", "mid", "post", "p_up", "p_down", "plambda"]
        est_fields = ["est_"+ f for f in fields]

        logger.info("(truth,estimated)")
        for field, est_field in zip(fields, est_fields):
            logger.info(f"{field}: {result[field], result[est_field]}")

        logger.info(f"AIC: {result['AIC']:.2f}")
        logger.info(f"ev_string: {result['ev_str']}")
        logger.info(f"est_ev_string: {result['est_ev_str']}")

    def process_results(self) -> List[Dict[str, Any]]:
        """
        Load the simulation results, add the test case number to them, print their info, and return them.

        Returns:
            A list of the results of the simulation, or an empty list if the file couldn't be loaded.
        """
        self._add_test_case_to_results()

        if self.data and "solutions" in self.data:
            logger.info(self.filename)
            logger.info(f"length of results: {len(self.data['solutions'])}")
            for result in self.data["solutions"]:
                self._print_result_info(result)
        else:
            return []

        return copy.deepcopy(self.data["solutions"])

def get_filenames_sorted_by_time(simulation_filename: str) -> List[str]:
    """
    Fetch all filenames associated with a given simulation and return them sorted by their creation times.

    Args:
        simulation_filename: The name of the simulation.

    Returns:
        A list of filenames sorted by their creation times.
    """
    file_pattern = f"SIMULATIONS/{simulation_filename}_*smaller.pickle"
    files = glob.glob(file_pattern)
    files.sort(key=os.path.getctime)
    return files

def ensure_sorted_by_worst_aic(results: List[Dict[str, Any]]) -> bool:
    """
    Check if the given list of simulation results is sorted by the 'AIC' field.

    Args:
        results: A list of simulation results.

    Returns:
        True if the results are sorted by 'AIC', False otherwise.
    """
    return all(results[i]['AIC'] >= results[i + 1]['AIC'] for i in range(len(results) - 1))

def process_files(filenames: List[str]) -> List[List[Dict[str, Any]]]:
    """
    Process a list of simulation result files and return their results.

    Args:
        filenames: A list of filenames.

    Returns:
        A list of lists, each containing the results of a simulation.
    """
    list_of_lists = []
    for filename in filenames:
        result = SimulationResult(filename).process_results()
        if result and ensure_sorted_by_worst_aic(result):
            list_of_lists.append(result)
    return list_of_lists

def get_all_best_estimates(list_of_lists: List[List[Dict[str, Any]]]) -> List[Dict[str, Any]]:
    """
    Get the best estimate from each list of simulation results.

    Args:
        list_of_lists: A list of lists of simulation results.

    Returns:
        A list of the best estimates.
    """
    return [x[-1] for x in list_of_lists if x]

def main(simulation_filename: str) -> None:
    """
    Main function of the script. Responsible for calling other helper functions and orchestrating
    the entire code flow.
    """
    sorted_filenames = get_filenames_sorted_by_time(simulation_filename)
    list_of_test_cases = process_files(sorted_filenames)
    best_estimates = get_all_best_estimates(list_of_test_cases)

    logger.info("BEST ESTIMATED")
    for estimate in best_estimates:
        logger.info(estimate)

