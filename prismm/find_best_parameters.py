import os
import pandas as pd
import shelve

BASE_FILE_FOLDER = "/vast/scratch/users/lmcintosh/CN_SV_SNV_evolutionary_likelihood/"
SIMULATIONS_FILE_FOLDER = BASE_FILE_FOLDER + "SIMULATIONS"
print(SIMULATIONS_FILE_FOLDER)




def get_simulation_files():
    simulation_files = [f for f in os.listdir(SIMULATIONS_FILE_FOLDER) if f.startswith("find_best_parameter_simulation_") and f.endswith(".txt.bak")]
    print(f"Simulation files: {simulation_files}")
    return simulation_files


def generate_path(pre, mid, post):
    path = ""

    if pre == -1 and (mid != -1 or post != -1):
        raise ValueError("Invalid input: pre is -1 while mid or post is not -1")

    if mid == -1 and post != -1:
        raise ValueError("Invalid input: mid is -1 while post is not -1")

    if pre >= 0:
        path += str(pre)

    if mid >= 0:
        path += "G" + str(mid)

    if post >= 0:
        path += "G" + str(post)

    return path

def load_simulation_data(suffix):
    file_name = f'{SIMULATIONS_FILE_FOLDER}/{suffix.split(".txt")[0]}.txt'
    with shelve.open(file_name) as d:
        try:
            pre = d['pre']
            mid = d['mid']
            post = d['post']
            path = generate_path(pre, mid, post)
            data = {
                'path': path,
                'max_number_of_solutions': d['max_number_of_solutions'],
                'max_default_path_length': d['max_default_path_length'],
                'prob_dist_filter': d['prob_dist_filter'],
                'path_length_diff': d['path_length_diff'],
                'searchable_likelihoods': d['searchable_likelihoods'],
                'marginal_likelihoods': d['marginal_likelihoods'],
                'top_likelihoods': d['top_likelihoods'],
                'pre': pre,
                'mid': mid,
                'post': post
            }
        except KeyError as e:
            print(f"KeyError occurred in file: {file_name}")
            print(f"Missing key: {e}")
            raise
    return data


def get_path_in_likelihoods(simulation_data, likelihoods_key):
    return simulation_data['path'] in simulation_data[likelihoods_key]['path'].values


def get_likelihoods_rank(simulation_data, likelihoods_key):
    if get_path_in_likelihoods(simulation_data, likelihoods_key):
        return simulation_data[likelihoods_key]['likelihood'].rank(ascending=False)[simulation_data[likelihoods_key]['path'] == simulation_data['path']].iloc[0]
    else:
        return None


def process_simulation_file(file_name):
    print(f"Processing {file_name}")
    simulation_data = load_simulation_data(file_name)
    test_case = int(file_name.split(".txt")[0].split("_")[-1])
    result = {
        'test_case': test_case,
        'max_number_of_solutions': simulation_data['max_number_of_solutions'],
        'max_default_path_length': simulation_data['max_default_path_length'],
        'prob_dist_filter': simulation_data['prob_dist_filter'],
        'path_length_diff': simulation_data['path_length_diff'],
        'path_in_searchable_likelihoods': get_path_in_likelihoods(simulation_data, 'searchable_likelihoods'),
        'searchable_likelihoods_rank': get_likelihoods_rank(simulation_data, 'searchable_likelihoods'),
        'path_in_marginal_likelihoods': get_path_in_likelihoods(simulation_data, 'marginal_likelihoods'),
        'marginal_likelihoods_rank': get_likelihoods_rank(simulation_data, 'marginal_likelihoods'),
        'path_in_top_likelihoods': get_path_in_likelihoods(simulation_data, 'top_likelihoods'),
        'top_likelihoods_rank': get_likelihoods_rank(simulation_data, 'top_likelihoods'),
        'path': simulation_data['path'],
        'pre': simulation_data['pre'],
        'mid': simulation_data['mid'],
        'post': simulation_data['post']
    }
    print(f"Result: {result}")
    return result



def generate_summary_df():
    simulation_files = get_simulation_files()
    data = [process_simulation_file(file_name) for file_name in simulation_files]
    return pd.DataFrame(data)


def save_summary_to_csv(df, file_name="simulation_summary.csv"):
    df.to_csv(file_name, index=False)


if __name__ == "__main__":
    summary_df = generate_summary_df()
    print(f"Summary dataframe:\n{summary_df}")
    save_summary_to_csv(summary_df)


