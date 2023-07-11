import collections
import numpy as np
import scipy.stats as stats
import typing as tp
import logging
import pandas as pd

# Set up logging.
logging.basicConfig(level=logging.INFO)


def empirical_distribution(distances: tp.List[float]) -> tp.List[tp.Tuple[float, float]]:
    """Calculate empirical distribution.

    Args:
        distances: A list of distances.

    Returns:
        A sorted list of tuples where each tuple contains a unique element from distances and its frequency.
    """
    counter = collections.Counter(distances)
    length = len(distances)
    return [(i, count / length) for i, count in sorted(counter.items())]


def calculate_summary_statistics(distances: tp.List[float]) -> tp.Dict[str, tp.Optional[float]]:
    """Calculate summary statistics.

    Args:
        distances: A list of data points.

    Returns:
        A dictionary containing summary statistics.
    """
    if not distances:
        print("Warning: no data to calculate summary statistics.")
        return {
            'mean': None,
            'std_dev': None,
            'min': None,
            '25th_percentile': None,
            'median': None,
            '75th_percentile': None,
            'max': None,
        }
    
    summary_stats = {
        'mean': np.mean(distances),
        'std_dev': np.std(distances),
        'min': np.min(distances),
        '25th_percentile': np.percentile(distances, 25),
        'median': np.median(distances),
        '75th_percentile': np.percentile(distances, 75),
        'max': np.max(distances),
    }
    return summary_stats


def print_summary_statistics(data_pairs: tp.List[tp.Tuple[float, float]], title: str) -> None:
    """Print summary statistics.

    Args:
        data_pairs: A list of tuples, each containing two paired data points.
        title: Title for the summary statistics.
    """
    logging.info(title)
    summary_stats = calculate_summary_statistics(data_pairs)
    for stat_name, stat_values in summary_stats.items():
        logging.info(f'{stat_name.capitalize()}: {stat_values}')
    logging.info('\n')
    return summary_stats


def summarize_p_similarity(all_results: tp.List[tp.Dict[str, tp.Any]]) -> None:
    """Summarize p similarity.

    Args:
        all_results: A list of dictionaries, each containing result data.
    """
    p_up_dists = [result["p_up"] - result["est_p_up"]/100.0 for result in all_results if result is not None]
    p_down_dists = [result["p_down"] - result["est_p_down"]/100.0 for result in all_results if result is not None]
   
    results = {}
    results["p_up_summary"] = print_summary_statistics(p_up_dists, "p_up diffs")
    results["p_down_summary"] = print_summary_statistics(p_down_dists, "p_down diffs")
    return results



def summarize_plambda_similarity(all_results: tp.List[tp.Dict[str, tp.Any]], label: str) -> dict:
    """Summarize plambda similarity for zero edit distance.

    Args:
        all_results: A list of dictionaries, each containing result data.
        label: A label to differentiate between different calls.
    """
    plambda_dists = []
    for result in all_results:
        if 'A' not in result['ev_str']:
            plambda_dists.append(0)
        else:
            plambda_dists.append(1 - result["est_plambda"] / result["plambda"])
    
    results = {
        f"plambda_summary_{label}": print_summary_statistics(plambda_dists, f"plambda diffs for {label}")
    }
    
    return results


def aggregate_similarity_scores(list_of_dfs: tp.List[pd.DataFrame]) -> pd.DataFrame:
    """Aggregate similarity scores.

    Args:
        list_of_dfs: A list of dataframes each having a 'Similarity Score' column.

    Returns:
        An aggregated dataframe.
    """
    num_rows = len(list_of_dfs[0])
    similarity_scores = [[] for _ in range(num_rows)]

    for df in list_of_dfs:
        for i, score in enumerate(df['Similarity Score']):
            # Only append the score if it is a number
            if pd.notnull(score):
                similarity_scores[i].append(score)

    stats_dict = {
        'Min': lambda scores: min(scores) if scores else np.nan,
        'Max': lambda scores: max(scores) if scores else np.nan,
        'Mean': lambda scores: np.mean(scores) if scores else np.nan,
        'Median': lambda scores: np.median(scores) if scores else np.nan,
        'StdDev': lambda scores: np.std(scores, ddof=0) if scores else np.nan,
        #'Mode': lambda scores: stats.mode(scores, nan_policy='omit')[0][0] if scores and len(set(scores)) > 1 else (scores[0] if scores and len(scores) > 0 else np.nan),
        'Mode': lambda scores: stats.mode(scores, nan_policy='omit', keepdims=True)[0][0] if scores and len(set(scores)) > 1 else (scores[0] if scores and len(scores) > 0 else np.nan),
        'IQR': lambda scores: stats.iqr(scores, nan_policy='omit') if scores else np.nan,
    }

    agg_df = list_of_dfs[0].copy()
    agg_df = agg_df.drop('Similarity Score', axis=1)
    agg_df.rename(columns={
        'Copy Number Distance Function': 'CN dist',
        'Epoch Created Distance Function': 'EC dist',
        'Normalising Constant': 'Norm Const'
    }, inplace=True)

    for stat_name, stat_func in stats_dict.items():
        agg_df[stat_name] = [stat_func(scores) for scores in similarity_scores]

    # Round values to 3 decimal places
    agg_df = agg_df.round(3)

    return agg_df

