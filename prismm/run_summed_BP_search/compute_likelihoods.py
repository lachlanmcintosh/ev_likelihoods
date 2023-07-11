import pandas as pd
import logging
from utils.path_code_to import path_code_to_path_length, path_code_num_anueploidy_epochs

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
    
# Define helper functions
def ensure_integer_columns(df, columns):
    for column in columns:
        if column in df.columns:
            df[column] = df[column].astype(int)
    return df

def get_likelihood_threshold(df):
    df['likelihood'] = pd.to_numeric(df['likelihood'])
    most_likely_likelihood = df['likelihood'].max()
    return most_likely_likelihood

def filter_by_likelihood(df, threshold):
    df = df[df['likelihood'] >= threshold]
    return df

def filter_by_path_length(df, x):
    df['path_length'] = df['path'].apply(path_code_to_path_length)
    df = df[df['path_length'] <= df['path_length'].min() + x]
    df.drop(columns=['path_length'], inplace=True)
    return df

def get_dataframe_subset(df, max_number_of_solutions, default_paths):
    df.sort_values('likelihood', ascending=False, inplace=True)
    df.drop_duplicates(subset='path', keep='first', inplace=True)
    result = df.head(max_number_of_solutions)
    
    for path in default_paths:
        path_rows = df[df["path"] == path]
        if not path_rows.empty:
            best_row = path_rows.loc[path_rows['likelihood'].idxmax()].to_frame().T
            result = pd.concat([result, best_row])
    
    result.drop_duplicates(subset='path', keep='first', inplace=True)
    result = result[result["likelihood"] > 0]
    
    return result

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


def create_searchable_likelihoods(marginal_likelihoods, top_likelihoods):
    marginal_likelihoods = marginal_likelihoods.reindex(columns=top_likelihoods.columns)
    searchable_likelihoods = pd.concat([top_likelihoods, marginal_likelihoods], ignore_index=True)
    searchable_likelihoods = searchable_likelihoods.drop_duplicates(subset=['path', 'p_up', 'p_down'], keep='first')
    searchable_likelihoods.sort_values(by='likelihood', ascending=False, inplace=True)

    logging.debug("The searchable likelihoods are:")
    logging.debug(searchable_likelihoods)
    logging.debug("paths to search: " + str(len(searchable_likelihoods.index)))

    return searchable_likelihoods


def compute_likelihoods(likelihoods, max_number_of_solutions, default_paths, prob_dist_filter, path_length_diff):
    likelihoods = ensure_integer_columns(likelihoods, ['p_up', 'p_down'])

    marginal_likelihoods = likelihoods_to_marginal_likelihoods(likelihoods, max_number_of_solutions, default_paths)
    marginal_likelihoods = ensure_integer_columns(marginal_likelihoods, ['p_up', 'p_down'])
    marginal_likelihoods = filter_by_likelihood(marginal_likelihoods, get_likelihood_threshold(marginal_likelihoods) * prob_dist_filter)
    marginal_likelihoods = filter_by_path_length(marginal_likelihoods, path_length_diff)

    top_likelihoods = likelihoods_to_best_likelihood_by_path(likelihoods, max_number_of_solutions, default_paths)
    top_likelihoods = ensure_integer_columns(top_likelihoods, ['p_up', 'p_down'])
    top_likelihoods = filter_by_likelihood(top_likelihoods, get_likelihood_threshold(top_likelihoods) * prob_dist_filter)
    top_likelihoods = filter_by_path_length(top_likelihoods, path_length_diff)

    searchable_likelihoods = create_searchable_likelihoods(marginal_likelihoods, top_likelihoods)
    searchable_likelihoods = ensure_integer_columns(searchable_likelihoods, ['p_up', 'p_down'])
    searchable_likelihoods = filter_out_p_up_p_down_equals_zero_row(searchable_likelihoods)

    return marginal_likelihoods, top_likelihoods, searchable_likelihoods
