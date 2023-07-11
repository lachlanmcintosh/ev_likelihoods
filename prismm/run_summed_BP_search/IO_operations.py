from utils.FILES import SIMULATIONS_FILE_FOLDER
import pickle
import os

def load_results_from_file(test_case, simulation_name):
    file_name = f'{SIMULATIONS_FILE_FOLDER}/{simulation_name}_{test_case}.pickle'
    with open(file_name, 'rb') as f:
        data = pickle.load(f)
    return data


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