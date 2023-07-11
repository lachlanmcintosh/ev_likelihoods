import glob
import os
import logging
from typing import List, Dict, Any, Optional
from dataclasses import dataclass

from clonal_trees.meta_analysis.simulation_result import SimulationResult


@dataclass
class SimulationProcessor:
    """
    A class for processing simulation results.

    Attributes:
        simulation_filename (str): The name of the simulation file to process.
    """
    simulation_filename: str

    def get_filenames_sorted_by_time(self) -> List[str]:
        """
        Retrieves the names of all simulation files sorted by their modification time.

        Returns:
            List[str]: A list of filenames sorted by modification time.
        """
        file_pattern = f"SIMULATIONS/{self.simulation_filename}_*smaller.pickle"
        files = glob.glob(file_pattern)
        files.sort(key=os.path.getctime)
        return files

    @staticmethod
    def ensure_sorted_by_worst_aic(results: List[Dict[str, Any]]) -> bool:
        """
        Checks if the results are sorted in descending order based on the 'AIC' key.

        Args:
            results (List[Dict[str, Any]]): A list of result dictionaries.

        Returns:
            bool: True if results are sorted by 'AIC' in descending order, False otherwise.

        Raises:
            ValueError: If the input list is empty.
        """
        if not results:
            raise ValueError("No results to sort.")
            
        return all(results[i]['AIC'] >= results[i + 1]['AIC'] for i in range(len(results) - 1))

    def process_files(self, filenames: List[str]) -> List[Optional[Dict[str, Any]]]:
        """
        Processes simulation files and ensures they are sorted in descending order based on the 'AIC' key.

        Args:
            filenames (List[str]): A list of filenames to process.

        Returns:
            List[Optional[Dict[str, Any]]]: A list of dictionaries containing processed results.
        """
        processed_results = []
        for filename in filenames:
            result = SimulationResult(filename).process_results()
            if result and self.ensure_sorted_by_worst_aic(result):
                processed_results.append(result)
        return processed_results

    @staticmethod
    def get_all_best_estimates(results: List[List[Dict[str, Any]]]) -> List[Dict[str, Any]]:
        """
        Extracts the best estimate from each list of results. 
        The best estimate is assumed to be the last element of each sorted list.

        Args:
            results (List[List[Dict[str, Any]]]): A list of lists of result dictionaries.

        Returns:
            List[Dict[str, Any]]: A list of dictionaries containing the best estimates.
        """
        return [result[-1] for result in results if result]

    def get_best_estimates(self) -> None:
        """
        Computes and logs the best estimates from all simulation files.
        """
        sorted_filenames = self.get_filenames_sorted_by_time()
        processed_files = self.process_files(sorted_filenames)
        best_estimates = self.get_all_best_estimates(processed_files)

        logging.info("BEST ESTIMATES")
        for estimate in best_estimates:
            logging.info(estimate)

