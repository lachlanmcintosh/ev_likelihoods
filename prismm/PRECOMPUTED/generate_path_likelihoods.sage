"""
This module precomputes matrices related to aneuploidy pathways based on user-provided parameters. 
The precomputed matrices are then used to calculate the likelihood of different pathways. 

Usage:
    python generate_path_likelihoods.py [p_up] [p_down] [max_CN] [path_length] [path_description]

Arguments:
    p_up: The probability of an upward mutation.
    p_down: The probability of a downward mutation.
    max_CN: The maximum copy number.
    path_length: The length of the path.
    path_description: The description of the path.
"""

import os
import argparse
import pickle
import numpy as np
from numpy import linalg as LA
import logging
from typing import List, Dict, Tuple, Optional

# Constants
MATRIX_DIRECTORY = "MATRICES"
GD_FILE = 'GD.sobj'
LOGGING_FORMAT = "%(levelname)s: %(message)s"

# Configure logging
logging.basicConfig(format=LOGGING_FORMAT, level=logging.INFO)


class MatrixPrecomputer:
    def __init__(self, 
                 p_up: float, 
                 p_down: float, 
                 max_CN: int, 
                 path_length: int, 
                 path_description: str):
        """
        Initialize the MatrixPrecomputer with the given parameters.

        Args:
            p_up (float): The probability of an upward mutation.
            p_down (float): The probability of a downward mutation.
            max_CN (int): The maximum copy number.
            path_length (int): The length of the path.
            path_description (str): The description of the path.
        """
        self.p_up = p_up
        self.p_down = p_down
        self.max_CN = max_CN
        self.path_length = path_length
        self.path_description = path_description

        # Output and intermediate filenames
        self.base_filename = f"{MATRIX_DIRECTORY}/subbed_mat_u{p_up}_d{p_down}_{path_description}.pickle"
        self.collated_output_file = f"{MATRIX_DIRECTORY}/collated_u{p_up}_d{p_down}_{path_description}.csv"
        self.powers_filename = self.base_filename.replace(".pickle", ".powers.pickle")
        self.precomputed_paths_filename = self.base_filename.replace(".pickle", ".precomputed_paths.pickle")

    def load_or_create_matrices(self) -> None:
        """
        Load or create matrices based on the provided parameters.

        Returns:
            None
        """
        if not os.path.isfile(self.collated_output_file):
            if not os.path.isfile(self.base_filename):
                self.create_substituted_matrix()
            self.anueploidy_matrix = self.load_matrix(self.base_filename)
            self.genome_doubling_matrix = self.load_matrix(GD_FILE)
            self.anueploidy_matrix, self.genome_doubling_matrix = self.ensure_correct_matrix_size()
            self.normalize_rows_in_anueploidy_matrix()
            self.all_paths = self.load_all_path_combinations()
            self.single_paths = self.separate_single_paths_from_all_paths()
            self.powers = self.calculate_powers_of_anueploidy_matrix_for_each_single_path()
            self.path_dict = self.precompute_paths()
            self.collate_path_likelihoods()


    def create_substituted_matrix(self):
        anueploidy_matrix = self.load_matrix(f"{MATRIX_DIRECTORY}/matrix_{self.path_description}.sobj")
        substituted_matrix = self.substitute_matrix(anueploidy_matrix)
        with open(self.base_filename, 'wb') as m_output:
            pickle.dump(substituted_matrix, m_output)

    def substitute_matrix(self, anueploidy_matrix: np.array) -> np.array:
        substituted_matrix = anueploidy_matrix.subs(u=self.p_up/100,d=self.p_down/100)
        substituted_matrix = substituted_matrix.apply_map(np.float64)
        substituted_matrix = substituted_matrix.numpy(dtype='double')
        return substituted_matrix

    def load_matrix(self, filename: str) -> np.array:
        with open(filename,'rb') as m_data:
            matrix = pickle.load(m_data)
        return matrix

    def ensure_correct_matrix_size(self) -> Tuple[np.array, np.array]:
        max_CN_plus_two = self.max_CN + 2
        anueploidy_matrix = self.anueploidy_matrix[:max_CN_plus_two,:max_CN_plus_two]
        genome_doubling_matrix = self.genome_doubling_matrix[:max_CN_plus_two,:max_CN_plus_two]
        return anueploidy_matrix, genome_doubling_matrix

    def normalize_rows_in_anueploidy_matrix(self):
        for row in range(self.max_CN + 1):
            total = np.sum(self.anueploidy_matrix[row,:])
            if total != 0:
                self.anueploidy_matrix[row,:] /= total

    def load_all_path_combinations(self) -> np.array:
        return self.load_matrix(f"all_path_combinations_{self.path_description}.sobj")

    def separate_single_paths_from_all_paths(self) -> List[str]:
        return [x for x in self.all_paths if "G" not in x]

    def calculate_powers_of_anueploidy_matrix_for_each_single_path(self) -> Dict[int, np.array]:
        if os.path.isfile(self.powers_filename):
            with open(self.powers_filename,'rb') as infile:
                powers = pickle.load(infile)
        else:
            powers = self.calculate_and_store_powers()
        return powers

    def calculate_and_store_powers(self) -> Dict[int, np.array]:
        powers = {}
        for path in self.single_paths:
            if int(path) not in powers:
                path_likelihood_mat = LA.matrix_power(self.anueploidy_matrix,int(path))
                powers[int(path)] = path_likelihood_mat
        with open(self.powers_filename,'wb') as infile:
            pickle.dump(powers,infile)
        return powers

    def precompute_paths(self) -> Dict[str, np.array]:
        if os.path.isfile(self.precomputed_paths_filename):
            try:
                with open(self.precomputed_paths_filename,'rb') as precomputed_data:
                    path_dict = pickle.load(precomputed_data)
            except Exception as e:
                logging.error(f"Failed to load precomputed paths due to: {str(e)}")
                os.remove(self.precomputed_paths_filename)
                path_dict = {}
        else:
            path_dict = self.precompute_paths_for_all_paths({})
        return path_dict

    def precompute_paths_for_all_paths(self, path_dict: Dict[str, np.array]) -> Dict[str, np.array]:
        count = 0
        for path in self.all_paths:
            if path in path_dict:
                continue
            path_dict[path] = self.precompute_path(path)
            if count % 200 == 0:
                with open(self.precomputed_paths_filename,'wb') as precomputed_data:
                    pickle.dump(path_dict, precomputed_data)
            count = count + 1
        return path_dict

    def precompute_path(self, path: str) -> np.array:
        splits = path.split("G")
        G1 = 0
        G2 = 0
        if len(splits) == 1:
            pre = int(path)
            mid = 0
            post = 0
        if len(splits) == 2:
            pre = int(splits[0])
            mid = 0
            post = int(splits[1])
            G1 = 1
        if len(splits) == 3:
            pre = int(splits[0])
            mid = int(splits[1])
            post = int(splits[2])
            G1 = 1
            G2 = 1

        # Compute the probabilities for this path
        path_likelihood_mat = self.powers[pre]

        if G1 > 0:
            path_likelihood_mat = np.matmul(path_likelihood_mat, self.genome_doubling_matrix)

        path_likelihood_mat = np.matmul(path_likelihood_mat, self.powers[mid])

        if G2 > 0:
            path_likelihood_mat = np.matmul(path_likelihood_mat, self.genome_doubling_matrix)

        path_likelihood_mat = np.matmul(path_likelihood_mat, self.powers[post])

        return path_likelihood_mat

    def collate_path_likelihoods(self):
        with open(self.precomputed_paths_filename, 'rb') as myd_data:
            path_dict = pickle.load(myd_data)
        out_text = self.collate_path_likelihoods_for_all_paths(path_dict)
        with open(self.collated_output_file,'w') as output_file:
            output_file.write(out_text)

    def collate_path_likelihoods_for_all_paths(self, path_dict: Dict[str, np.array]) -> str:
        out_text = ""
        for path in path_dict:
            path_likelihood_mat = path_dict[path]
            line = f"{self.p_up},{self.p_down},{path},"
            line += ",".join([str(x) for x in path_likelihood_mat[1,:]]) + "\n"
            out_text += line
        return out_text

def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments.

    Returns:
        argparse.Namespace: Namespace with parsed arguments.
    """
    parser = argparse.ArgumentParser(description="Precompute matrices for anueploidy pathways.")
    parser.add_argument('p_up', type=float, help='The probability of an upward mutation.')
    parser.add_argument('p_down', type=float, help='The probability of a downward mutation.')
    parser.add_argument('max_CN', type=int, help='The maximum copy number.')
    parser.add_argument('path_length', type=int, help='The length of the path.')
    parser.add_argument('path_description', type=str, help='The description of the path.')
    return parser.parse_args()


def main() -> None:
    """
    The main function that initializes and runs the MatrixPrecomputer.

    Returns:
        None
    """
    args = parse_args()
    MatrixPrecomputer(args.p_up, args.p_down, args.max_CN, args.path_length, args.path_description).load_or_create_matrices()


if __name__ == "__main__":
    main()

