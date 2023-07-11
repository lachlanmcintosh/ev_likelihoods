import pandas as pd
import logging


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
