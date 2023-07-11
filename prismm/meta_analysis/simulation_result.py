import os
import glob
import pickle
import re
import copy
import logging
from typing import List, Optional, Any, Dict, Union

from fields import Field

class SimulationResult:
    """
    A class to represent a Simulation Result.

    Attributes
    ----------
    filename : str
        the name of the file to load data from
    data : dict
        the data loaded from the file

    Methods
    -------
    _load_file():
        Loads a pickle file and returns a dictionary. Handles error during file loading.
    _get_test_case():
        Returns an integer representing the test case parsed from the filename. If no test case is found, raise ValueError.
    _add_test_case_to_results():
        Adds the test case to the result data.
    _print_result_info(result: Dict[str, Any]):
        Logs the information about the result.
    process_results():
        Processes the result data and logs their information. Returns a deep copy of the solutions in the data.
    """

    def __init__(self, filename: str):
        """
        Constructs all the necessary attributes for the SimulationResult object.
        """
        self.filename = filename
        self.data = self._load_file()

    def _load_file(self) -> Optional[dict]:
        """
        Loads a pickle file and returns a dictionary.
        If an error occurs during the process, log the error and return None.
        """
        try:
            with open(self.filename, 'rb') as file:
                return pickle.load(file)
        except EOFError:
            logging.warning(f"EOFError: {self.filename} is incomplete. The file might still be generating.")
        except Exception as e:
            logging.error(f"An error occurred while loading {self.filename}: {e}")
        return None

    def _get_test_case(self) -> int:
        """
        Parses the test case number from the filename.
        If no test case is found, raise a ValueError.
        """
        match = re.search(r'(\d+)', self.filename)
        if match is None:
            raise ValueError(f"Test case not found in filename: {self.filename}")
        return int(match.group())

    def _add_test_case_to_results(self) -> None:
        """
        Adds the parsed test case to each solution in the result data.
        """
        test_case = self._get_test_case()

        if self.data and "solutions" in self.data:
            for result in self.data["solutions"]:
                result['test_case'] = test_case

    def _print_result_info(self, result: Dict[str, Any]) -> None:
        """
        Logs the information about the result including truth, estimated values, AIC, and EV strings.
        """
        est_fields = ["est_" + f.value for f in Field]

        logging.info("(truth,estimated)")
        for field, est_field in zip(list(Field), est_fields):
            logging.info(f"{field.value}: {result[field.value], result[est_field]}")

        logging.info(f"AIC: {result['AIC']:.2f}")
        logging.info("ev_string: " + str(result["ev_str"]))
        logging.info("est_ev_string: " + str(result["est_ev_str"]))

    def process_results(self) -> List[Dict[str, Union[str, float]]]:
        """
        Processes the result data by adding test case information and logging their details.
        Returns a deep copy of the solutions in the data if they exist, else returns an empty list.
        """
        self._add_test_case_to_results()

        if self.data and "solutions" in self.data:
            logging.info(self.filename)
            logging.info(f"length of results: {len(self.data['solutions'])}")
            for result in self.data["solutions"]:
                self._print_result_info(result)
            return copy.deepcopy(self.data["solutions"])
        else:
            return []

