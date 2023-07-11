import numpy as np
import sys
import re
import pandas as pd
from typing import Dict
from constants import *
from general_functions import *
import logging
import pickle
import os
import argparse



def calculate_likelihoods(all_data, observed_CN_multiplicities):
    observed_CN_multiplicities_str = {
        str(key): value for key, value in observed_CN_multiplicities.items()
    }
    lls = np.sum(
        [all_data[copy] * multiplicity for copy, multiplicity in observed_CN_multiplicities_str.items()],
        axis=0
    )
    
    assert np.shape(lls) == np.shape(all_data["0"])
    likelihoods = np.exp(lls)
    if "0" in observed_CN_multiplicities:
        likelihoods = likelihoods / (1 - np.exp(all_data["0"]))**observed_CN_multiplicities["0"]
    return likelihoods


def test_calculate_likelihoods():
    # Test Case 1
    all_data = {
        "0": np.array([0.0, 0.0, 0.0]),
        "1": np.array([0.5, 0.5, 0.5]),
        "2": np.array([1.0, 1.0, 1.0]),
        "3": np.array([1.5, 1.5, 1.5])
    }
    observed_CN_multiplicities = {
        0: 2,
        1: 1,
        2: 3
    }
    
    expected_likelihoods = np.array([0.01831564, 0.01831564, 0.01831564])
    
    result = calculate_likelihoods(all_data, observed_CN_multiplicities)
    assert np.allclose(result, expected_likelihoods), f"Expected {expected_likelihoods}, but got {result}"

    # Test Case 2
    all_data = {
        "0": np.array([0.1, 0.2, 0.3]),
        "1": np.array([0.4, 0.5, 0.6]),
        "2": np.array([0.7, 0.8, 0.9]),
        "3": np.array([1.0, 1.1, 1.2])
    }
    observed_CN_multiplicities = {
        1: 2,
        2: 1
    }
    
    expected_likelihoods = np.array([1.778271, 2.178106, 2.576879])
    
    result = calculate_likelihoods(all_data, observed_CN_multiplicities)
    assert np.allclose(result, expected_likelihoods, rtol=1e-5), f"Expected {expected_likelihoods}, but got {result}"

    # Test Case 3
    all_data = {
        "0": np.array([0.2, 0.3, 0.4]),
        "1": np.array([0.5, 0.6, 0.7]),
        "2": np.array([0.8, 0.9, 1.0]),
        "3": np.array([1.1, 1.2, 1.3])
    }
    observed_CN_multiplicities = {
        0: 1,
        1: 1,
        2: 1,
        3: 1
    }
    
    expected_likelihoods = np.array([1.513225, 1.878165, 2.218444])
    
    result = calculate_likelihoods(all_data, observed_CN_multiplicities)
    assert np.allclose(result, expected_likelihoods, rtol=1e-5), f"Expected {expected_likelihoods}, but got {result}"


#test_calculate_likelihoods()


def read_pickle_with_custom_columns(file_name):
    all_data = pd.read_pickle(file_name)
    return all_data

def CN_multiplicities_to_likelihoods(observed_CN_multiplicities: Dict[int, int]):
    file_name = PRECOMPUTED_FILE_FOLDER + "collated_p8_v4_logged.pkl"
    p_value = int(re.findall(r"_p(\d+)_", file_name)[0])
    all_data = read_pickle_with_custom_columns(file_name)
    likelihoods = calculate_likelihoods(all_data, observed_CN_multiplicities)
    named_likelihoods = all_data[["p_up", "p_down", "path"]].copy()
    named_likelihoods.insert(3, "likelihood", likelihoods, True)
    named_likelihoods.replace([np.inf, -np.inf], np.nan, inplace=True)
    named_likelihoods.dropna(axis=0, inplace=True)
    named_likelihoods.sort_values(by=['likelihood'], inplace=True, ascending=False)
    total = np.nansum(named_likelihoods["likelihood"])
    named_likelihoods["likelihood"] /= total
    return named_likelihoods

def test_CN_multiplicities_to_likelihoods():
    def read_pickle_with_custom_columns(file_name):
        # Dummy implementation
        return pd.DataFrame({
            "0": [0, 1, 2],
            "1": [0, 1, 2],
            "2": [0, 1, 2],
            "3": [0, 1, 2],
            "p_up": [0.1, 0.2, 0.3],
            "p_down": [0.3, 0.2, 0.1],
            "path": ["path1", "path2", "path3"],
        })

    observed_CN_multiplicities = {
        0: 2,
        1: 1,
        2: 3
    }

    expected_named_likelihoods = pd.DataFrame({
        "p_up": [0.1, 0.2, 0.3],
        "p_down": [0.3, 0.2, 0.1],
        "path": ["path1", "path2", "path3"],
        "likelihood": [0.31830989, 0.31830989, 0.36338023]
    })

    result = CN_multiplicities_to_likelihoods(observed_CN_multiplicities)
    pd.testing.assert_frame_equal(result, expected_named_likelihoods, check_exact=False, rtol=1e-5)

    # Additional test cases can be added here with different observed_CN_multiplicities values and read_pickle_with_custom_columns function

#test_CN_multiplicities_to_likelihoods()


def likelihoods_to_marginal_likelihoods(likelihoods, max_number_of_solutions, default_paths):
    marginal_likelihoods = likelihoods.copy()

    marginal_likelihoods["mean_p_up"] = marginal_likelihoods["p_up"] * marginal_likelihoods["likelihood"]
    marginal_likelihoods["mean_p_down"] = marginal_likelihoods["p_down"] * marginal_likelihoods["likelihood"]

    marginal_likelihoods = marginal_likelihoods.groupby(["path"], as_index=False)[['likelihood', 'mean_p_up', 'mean_p_down']].sum()

    marginal_likelihoods.sort_values(by=['likelihood'], inplace=True, ascending=False)
    result = marginal_likelihoods.head(max_number_of_solutions)

    for path in default_paths:
        best_row = marginal_likelihoods.loc[marginal_likelihoods["path"] == path]
        result = pd.concat([result, best_row])

    result = result.drop_duplicates(subset='path', keep='first')

    result["mean_p_up"] /= result["likelihood"]
    result["mean_p_down"] /= result["likelihood"]

    result.sort_values(by=['likelihood'], inplace=True, ascending=False)

    result["p_up"] = round(result["mean_p_up"])
    result["p_down"] = round(result["mean_p_down"])

    result = result[["p_up", "p_down", "path", "likelihood"]]

    # Remove rows with a likelihood of 0
    result = result.loc[result["likelihood"] > 0]

    return result



def test_likelihoods_to_marginal_likelihoods():
    likelihoods = pd.DataFrame({
        "p_up": [0.1, 0.2, 0.3, 0.4],
        "p_down": [0.3, 0.2, 0.1, 0.0],
        "path": ["path1", "path2", "path3", "path4"],
        "likelihood": [0.25, 0.25, 0.25, 0.25]
    })

    max_number_of_solutions = 2
    default_paths = ["path1", "path4"]

    expected_result = pd.DataFrame({
        "p_up": [0.0, 0.0, 1.0, 4.0],
        "p_down": [1.0, 0.0, 3.0, 0.0],
        "path": ["path4", "path1", "path2", "path3"],
        "likelihood": [0.25, 0.25, 0.25, 0.25]
    })

    result = likelihoods_to_marginal_likelihoods(likelihoods, max_number_of_solutions, default_paths)
    pd.testing.assert_frame_equal(result, expected_result)

    # Additional test cases can be added here with different likelihoods, max_number_of_solutions, and default_paths values


#test_likelihoods_to_marginal_likelihoods()



def likelihoods_to_best_likelihood_by_path(likelihoods, max_number_of_solutions, default_paths):
    df = likelihoods.copy()
    
    df = df.sort_values('likelihood', ascending=False)
    
    df = df.drop_duplicates(subset='path', keep='first')
    
    result = df.head(max_number_of_solutions)
    
    for path in default_paths:
        path_rows = df.loc[df["path"] == path]
        if path_rows.empty:
            continue
        best_row = path_rows.loc[path_rows['likelihood'].idxmax()].to_frame().T
        result = pd.concat([result, best_row])
    
    result = result.drop_duplicates(subset='path', keep='first')
    
    result.sort_values('likelihood', ascending=False, inplace=True)
    result = result.loc[result["likelihood"] > 0]
    
    return result


def test_likelihoods_to_best_likelihood_by_path():
    likelihoods = pd.DataFrame({
        "p_up": [0.1, 0.2, 0.3, 0.4],
        "p_down": [0.3, 0.2, 0.1, 0.0],
        "path": ["path1", "path2", "path3", "path4"],
        "likelihood": [0.1, 0.2, 0.3, 0.4]
    })

    max_number_of_solutions = 2
    default_paths = ["path1", "path4"]

    expected_result = pd.DataFrame({
        "p_up": [0.4, 0.3, 0.1, 0.2],
        "p_down": [0.0, 0.1, 0.3, 0.2],
        "path": ["path4", "path3", "path1", "path2"],
        "likelihood": [0.4, 0.3, 0.1, 0.2]
    })

    result = likelihoods_to_best_likelihood_by_path(likelihoods, max_number_of_solutions, default_paths)
    pd.testing.assert_frame_equal(result, expected_result)

    # Additional test cases can be added here with different likelihoods, max_number_of_solutions, and default_paths values


# test_likelihoods_to_best_likelihood_by_path()


def load_results_from_file(test_case, simulation_name):
    file_name = f'{SIMULATIONS_FILE_FOLDER}/{simulation_name}_{test_case}.pickle'
    with open(file_name, 'rb') as f:
        data = pickle.load(f)
    return data

def generate_default_paths(max_default_path_length):
    default_paths = [str(x) for x in range(max_default_path_length)]
    default_paths += [str(x) + "G" + str(y)
                      for x in range(max_default_path_length)
                      for y in range(max_default_path_length)
                      if x + y <= max_default_path_length]

    return default_paths



def print_path_likelihoods(likelihoods, searchable_likelihoods, marginal_likelihoods, top_likelihoods, default_paths, data):
    def print_dataframe_rows(df, label, path):
        path_in_df = any(df['path'] == path)
        if path_in_df:
            path_rows = df.loc[df['path'] == path]
            logging.info(f"{label} for path {path}:\n{path_rows}")
        else:
            logging.info(f"Path {path} not found in {label}")
            
    # Create the path using the data dictionary values for pre, mid, and post
    path = generate_path(data['pre'], data['mid'], data['post'])
    
    # Print the pre, mid, and post values and the generated path
    logging.info(f"pre: {data['pre']}, mid: {data['mid']}, post: {data['post']} => path: {path}")
    
    # Print the likelihood row for that path
    path_in_likelihoods = any(likelihoods['path'] == path)
    if path_in_likelihoods:
        print_dataframe_rows(likelihoods, "Likelihood", path)
    else:
        raise ValueError(f"Path {path} not found in likelihoods")
        
    # Print the rows for the given path in other dataframes
    print_dataframe_rows(marginal_likelihoods, "Marginal likelihood", path)
    print_dataframe_rows(top_likelihoods, "Top likelihood", path)
    print_dataframe_rows(searchable_likelihoods, "Searchable likelihood", path)
    
    # Check if the path is in default_paths and print a message accordingly
    logging.info(f"Default paths searched through are {default_paths}")
    if path in default_paths:
        logging.info(f"Path {path} is a default path")
    else:
        logging.info(f"Path {path} not found in default paths")

    # Calculate the marginal likelihood for the path
    path_likelihood_rows = likelihoods.loc[likelihoods['path'] == path]
    marginal_likelihood = path_likelihood_rows['likelihood'].sum()
    logging.info(f"Marginal likelihood for path {path}: {marginal_likelihood}")

    # Calculate the number of unique paths with a higher likelihood and get the paths
    count_higher_likelihood_paths(likelihoods=marginal_likelihoods, path=path, name="marginal likelihoods")

    # same for top likelihoods
    count_higher_likelihood_paths(likelihoods=top_likelihoods, path=path, name="top likelihoods")


def count_higher_likelihood_paths1(likelihoods, path):
  """
  This function counts the number of paths in the `likelihoods` dataframe with higher likelihood than the path passed in as an argument. It also returns the unique paths with higher likelihood and the ratio of the likelihood of the most likely path to the likelihood of the path passed in as an argument.

  Args:
    likelihoods: A Pandas DataFrame containing the likelihoods of all paths.
    path: The path to count the number of higher likelihood paths for.

  Returns:
    A tuple containing the number of higher likelihood paths, the unique paths with higher likelihood, and the ratio of the likelihood of the most likely path to the likelihood of the path passed in as an argument.
  """

  # Get the row with the given path
  path_row = likelihoods.loc[likelihoods['path'] == path]

  # If the row does not exist, create it and set the likelihood to 0
  if path_row.empty:
    path_row = pd.DataFrame({'path': [path], 'likelihood': [0]})
    path_likelihood = 0
  else:
    path_likelihood = path_row['likelihood'].values[0]

  # Get the rows with higher likelihood
  higher_likelihood_rows = likelihoods[likelihoods['likelihood'] > path_likelihood]

  # Get the unique paths with higher likelihood
  unique_higher_likelihood_paths = higher_likelihood_rows['path'].unique()

  # Get the likelihood of the most likely path
  most_likely_likelihood = likelihoods['likelihood'].max()

  # Calculate the ratio of likelihoods
  ratio_of_likelihoods = most_likely_likelihood / path_likelihood

  return len(unique_higher_likelihood_paths), unique_higher_likelihood_paths, ratio_of_likelihoods

def count_higher_likelihood_paths2(likelihoods, path):
    # Get the row with the given path
    path_row = likelihoods.loc[likelihoods['path'] == path]

    # Check if path_row is empty
    if path_row.empty:
        return 0, [], float('nan')

    # Get the likelihood value for the given path
    path_likelihood = path_row['likelihood'].values[0]

    # Get the rows with higher likelihood
    higher_likelihood_rows = likelihoods[likelihoods['likelihood'] > path_likelihood]

    # Get the unique paths with higher likelihood
    unique_higher_likelihood_paths = higher_likelihood_rows['path'].unique()

    # Find the maximum likelihood among the higher likelihood paths
    max_likelihood = higher_likelihood_rows['likelihood'].max()

    # Calculate the likelihood ratio
    likelihood_ratio = max_likelihood / path_likelihood

    return len(unique_higher_likelihood_paths), unique_higher_likelihood_paths, likelihood_ratio


def count_higher_likelihood_paths(likelihoods, path, name):
    C1,P1,LR1 = count_higher_likelihood_paths1(likelihoods, path)
    C2,P2,LR2 = count_higher_likelihood_paths2(likelihoods, path)
    #assert(C1==C2), f"C1: {C1}, C2: {C2}"
    #assert(sorted(P1)==sorted(P2)), f"Sorted P1: {sorted(P1)}, Sorted P2: {sorted(P2)}"
    #assert(LR1==LR2), f"LR1: {LR1}, LR2: {LR2}"
    logging.info(f"Number of unique paths with higher {name} than path {path}: {C1}")
    logging.info(f"Unique paths with higher {name} than path {path}: {P1}")
    logging.info(f"Likelihood ratio of top path to {path}: {LR1}")
    return C1,P1,LR1


def pretty_print_data(data):
    for key, value in data.items():
        pretty_print(key)
        pretty_print(value)


def ensure_integer_columns(df, columns):
    for column in columns:
        if column in df.columns:
            df[column] = df[column].astype(int)
    return df

def compute_likelihoods(likelihoods, max_number_of_solutions, default_paths, prob_dist_filter, path_length_diff):
    likelihoods = ensure_integer_columns(likelihoods, ['p_up', 'p_down'])

    marginal_likelihoods = likelihoods_to_marginal_likelihoods(
        likelihoods=likelihoods,
        max_number_of_solutions=max_number_of_solutions,
        default_paths=default_paths
    )
    marginal_likelihoods = ensure_integer_columns(marginal_likelihoods, ['p_up', 'p_down'])

    top_likelihoods = likelihoods_to_best_likelihood_by_path(
        likelihoods=likelihoods,
        max_number_of_solutions=max_number_of_solutions,
        default_paths=default_paths
    )
    top_likelihoods = ensure_integer_columns(top_likelihoods, ['p_up', 'p_down'])

    # Filter marginal_likelihoods based on the likelihood values
    marginal_likelihood_threshold = get_likelihood_threshold(marginal_likelihoods)
    marginal_likelihoods = filter_rows_by_likelihood(marginal_likelihoods, marginal_likelihood_threshold * prob_dist_filter)

    # Filter top_likelihoods based on the likelihood values
    top_likelihood_threshold = get_likelihood_threshold(top_likelihoods)
    top_likelihoods = filter_rows_by_likelihood(top_likelihoods, top_likelihood_threshold * prob_dist_filter)

    # Apply filter_rows_by_path_length on marginal_likelihoods
    marginal_likelihoods = filter_rows_by_path_length(marginal_likelihoods, path_length_diff)

    # Apply filter_rows_by_path_length on top_likelihoods
    top_likelihoods = filter_rows_by_path_length(top_likelihoods, path_length_diff)

    searchable_likelihoods = create_searchable_likelihoods(marginal_likelihoods, top_likelihoods)
    searchable_likelihoods = ensure_integer_columns(searchable_likelihoods, ['p_up', 'p_down'])

    return marginal_likelihoods, top_likelihoods, searchable_likelihoods


def create_searchable_likelihoods(marginal_likelihoods, top_likelihoods):
    marginal_likelihoods = marginal_likelihoods.reindex(columns=top_likelihoods.columns)
    searchable_likelihoods = pd.concat([top_likelihoods, marginal_likelihoods], ignore_index=True)
    searchable_likelihoods = searchable_likelihoods.drop_duplicates(subset=['path', 'p_up', 'p_down'], keep='first')
    searchable_likelihoods.sort_values(by='likelihood', ascending=False, inplace=True)

    logging.debug("The searchable likelihoods are:")
    logging.debug(searchable_likelihoods)
    logging.debug("paths to search: " + str(len(searchable_likelihoods.index)))

    return searchable_likelihoods




def filter_out_p_up_p_down_equals_zero_row(df: pd.DataFrame) -> pd.DataFrame:
    df['aneuploidy_epochs'] = df['path'].apply(path_code_num_anueploidy_epochs)

    zero_p_up_down = df.loc[(df['p_up'] == 0) & (df['p_down'] == 0), 'aneuploidy_epochs']

    if zero_p_up_down.empty:
        min_idx = None
    else:
        min_idx = zero_p_up_down.idxmin()

    if min_idx is not None:
        min_row = df.loc[[min_idx]]
        df_filtered = df.loc[~((df['p_up'] == 0) & (df['p_down'] == 0) & (df.index != min_idx))]
    else:
        df_filtered = df.loc[~((df['p_up'] == 0) & (df['p_down'] == 0))]

    #df_filtered.drop(columns=['aneuploidy_epochs'], inplace=True)
    df_filtered = df_filtered.drop(columns=['aneuploidy_epochs'])


    return df_filtered


def test_filter_out_p_up_p_down_equals_zero_row():
    data = {
        "p_up": [33, 2, 97, 1, 1, 82, 0, 0, 0, 0],
        "p_down": [33, 2, 2, 1, 1, 13, 0, 0, 0, 0],
        "path": ["0G0", "1G0", "1", "0G1", "2G0", "0G0", "3G3", "3G4", "4G3", "5G0"],
        "likelihood": [0.989199, 0.001347, 0.001347, 0.000524, 0.000523, 0.000192, 0.000192, 0.000192, 0.000192, 0.000192]
    }

    df = pd.DataFrame(data)
    filtered_df = filter_out_p_up_p_down_equals_zero_row(df)

    expected_filtered_data = {
        "p_up": [33, 2, 97, 1, 1, 82, 0],
        "p_down": [33, 2, 2, 1, 1, 13, 0],
        "path": ["0G0", "1G0", "1", "0G1", "2G0", "0G0", "5G0"],
        "likelihood": [0.989199, 0.001347, 0.001347, 0.000524, 0.000523, 0.000192, 0.000192]
    }

    expected_filtered_df = pd.DataFrame(expected_filtered_data)

    # Reset the index of both DataFrames
    filtered_df.reset_index(drop=True, inplace=True)
    expected_filtered_df.reset_index(drop=True, inplace=True)

    # Compare the filtered DataFrame with the expected output
    try:
        pd.testing.assert_frame_equal(filtered_df, expected_filtered_df)
    except AssertionError as e:
        print("Expected:")
        print(expected_filtered_df)
        print("Actual:")
        print(filtered_df)
        raise e


test_filter_out_p_up_p_down_equals_zero_row()


def path_length(path: str) -> int:
    return len(path.replace("G", ""))



def filter_rows_by_path_length(df: pd.DataFrame, x: int) -> pd.DataFrame:
    df['path_length'] = df['path'].apply(path_code_to_path_length)

    # Convert the 'likelihood' column to a numeric dtype
    df['likelihood'] = pd.to_numeric(df['likelihood'])

    min_path_length = df['path_length'].min()
    most_likely_solution_path_length = df.loc[df['likelihood'].idxmax(), 'path_length']
    min_path_length = max(min_path_length, most_likely_solution_path_length)

    df_filtered = df.loc[df['path_length'] <= min_path_length + x]

    #df_filtered.drop(columns=['path_length'], inplace=True)
    df_filtered = df_filtered.drop(columns=['path_length'])

    return df_filtered

def test_filter_rows_by_path_length():
    data = {
        "p_up": [33, 2, 97, 1, 1, 82, 0, 0, 0, 0],
        "p_down": [33, 2, 2, 1, 1, 13, 0, 0, 0, 0],
        "path": ["0G0", "1G0", "1", "0G1", "2G0", "0G0", "3G3", "3G4", "4G3", "5G0"],
        "likelihood": [0.989199, 0.001347, 0.001347, 0.000524, 0.000523, 0.000192, 0.000192, 0.000192, 0.000192, 0.000192]
    }

    df = pd.DataFrame(data)
    x = 2
    filtered_df = filter_rows_by_path_length(df, x)

    expected_filtered_data = {
        "p_up": [33, 2, 97, 1, 1, 82, 0],
        "p_down": [33, 2, 2, 1, 1, 13, 0],
        "path": ["0G0", "1G0", "1", "0G1", "2G0", "0G0", "5G0"],
        "likelihood": [0.989199, 0.001347, 0.001347, 0.000524, 0.000523, 0.000192, 0.000192]
    }

    expected_filtered_df = pd.DataFrame(expected_filtered_data)

    # Reset the index of both DataFrames
    filtered_df.reset_index(drop=True, inplace=True)
    expected_filtered_df.reset_index(drop=True, inplace=True)

    # Compare the filtered DataFrame with the expected output
    pd.testing.assert_frame_equal(filtered_df, expected_filtered_df)

#test_filter_rows_by_path_length()


def filter_rows_by_likelihood(df: pd.DataFrame, likelihood_threshold: float) -> pd.DataFrame:
    df['likelihood'] = pd.to_numeric(df['likelihood'])
    df_filtered = df.loc[df['likelihood'] >= likelihood_threshold ]
    return df_filtered

def get_likelihood_threshold(df: pd.DataFrame) -> float:
    df['likelihood'] = pd.to_numeric(df['likelihood'])
    most_likely_likelihood = df['likelihood'].max()
    return most_likely_likelihood

def print_dataframes(dataframes: dict):
    for name, df in dataframes.items():
        logging.info(f"{name}:\n{df}\n")



def save_likelihoods(test_case, simulation_name, likelihoods, marginal_likelihoods, top_likelihoods, searchable_likelihoods, max_number_of_solutions, max_default_path_length, prob_dist_filter, path_length_diff):
    file_name = f'{SIMULATIONS_FILE_FOLDER}/{simulation_name}_{test_case}.pickle'
    
    data_to_dump = {
        'likelihoods': likelihoods,
        'marginal_likelihoods': marginal_likelihoods,
        'top_likelihoods': top_likelihoods,
        'searchable_likelihoods': searchable_likelihoods,
        'max_number_of_solutions': max_number_of_solutions,
        'max_default_path_length': max_default_path_length,
        'prob_dist_filter': prob_dist_filter,
        'path_length_diff': path_length_diff
    }
    
    if os.path.exists(file_name):
        with open(file_name, 'rb') as f:
            old_data = pickle.load(f)
        old_data.update(data_to_dump)
        data_to_dump = old_data
    
    with open(file_name, 'wb') as f:
        pickle.dump(data_to_dump, f)



def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Process and compute likelihoods for a given test case."
    )

    parser.add_argument(
        "-t",
        "--test_case",
        type=int,
        help="The name or number of the test case you want to process."
    )

    parser.add_argument(
        "-f",
        "--simulation_filename",
        required=True,
        help="Name of the simulation."
    )

    parser.add_argument(
        "-s",
        "--max_number_of_solutions",
        type=int,
        default=20,
        help="Maximum number of solutions to compute."
    )

    parser.add_argument(
        "-d",
        "--max_default_path_length",
        type=int,
        default=3,
        help="Maximum length of the default paths."
    )

    parser.add_argument(
        "-p",
        "--prob_dist_filter",
        type=float,
        default=0.1,
        help="Probability distribution filter value."
    )

    parser.add_argument(
        "-l",
        "--path_length_diff",
        type=int,
        default=2,
        help="Path length difference value."
    )

    parser.add_argument(
        "-b",
        "--debug",
        action='store_true',
        help="Enable debug logging."
    )

    return parser.parse_args()



def main():
    args = parse_arguments()

    if args.debug:
        logging.basicConfig(level = logging.DEBUG)
    else:
        logging.basicConfig(level = logging.INFO)

    default_paths = generate_default_paths(max_default_path_length = args.max_default_path_length)

    data = load_results_from_file(test_case = args.test_case, simulation_name = args.simulation_filename)

    #pretty_print_data(data = data)

    likelihoods = CN_multiplicities_to_likelihoods(observed_CN_multiplicities = data['observed_CN_multiplicities'])
    marginal_likelihoods, top_likelihoods, searchable_likelihoods = compute_likelihoods(
        likelihoods = likelihoods,
        max_number_of_solutions = args.max_number_of_solutions,
        default_paths = default_paths,
        prob_dist_filter = args.prob_dist_filter,
        path_length_diff = args.path_length_diff
    )

    searchable_likelihoods = filter_out_p_up_p_down_equals_zero_row(searchable_likelihoods)

    # Dictionary of DataFrames with their names
    dataframes = {
        "Marginal Likelihoods" : marginal_likelihoods,
        "Top Likelihoods" : top_likelihoods,
        "Likelihoods" : likelihoods,
        "Searchable Likelihoods" : searchable_likelihoods
    }
    print_dataframes(dataframes)
    logging.info("finished_print_dataframes")

    save_likelihoods(
        test_case = args.test_case,
        simulation_name = args.simulation_filename,
        likelihoods = likelihoods,
        marginal_likelihoods = marginal_likelihoods,
        top_likelihoods = top_likelihoods,
        searchable_likelihoods = searchable_likelihoods,
        max_number_of_solutions = args.max_number_of_solutions,
        max_default_path_length = args.max_default_path_length,
        prob_dist_filter = args.prob_dist_filter,
        path_length_diff = args.path_length_diff
    )

    print_path_likelihoods(
        likelihoods = likelihoods,
        searchable_likelihoods = searchable_likelihoods,
        marginal_likelihoods = marginal_likelihoods,
        top_likelihoods = top_likelihoods,
        default_paths = default_paths,
        data = data
    )



if __name__ == "__main__":
    main()

