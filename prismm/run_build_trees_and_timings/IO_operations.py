from utils.FILES import SIMULATIONS_FILE_FOLDER
import os
import pickle


def save_dict_to_file(dictionary, test_case, simulation_filename):
    filename = f'{SIMULATIONS_FILE_FOLDER}/{simulation_filename}_{test_case}.pickle'
    
    if os.path.exists(filename):
        with open(filename, 'rb') as f:
            old_data = pickle.load(f)
        old_data.update(dictionary)
        dictionary = old_data
    
    with open(filename, 'wb') as f:
        pickle.dump(dictionary, f)


def load_results_from_file(test_case, simulation_filename):
    file_name = f'{SIMULATIONS_FILE_FOLDER}/{simulation_filename}_{test_case}.pickle'
    with open(file_name, 'rb') as f:
        data = pickle.load(f)
    return data
