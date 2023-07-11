import pickle
import os


def load_simulation_data(test_case, simulation_filename):
    """
    This function opens a file containing simulation data and returns the simulation results.

    Args:
    test_case (str): A string specifying which test case to load.
    simulation_filename (str): A string specifying the name of the simulation.

    Returns:
    dict: A dictionary containing the results of the simulation.
    """
    with open(f'{SIMULATIONS_FILE_FOLDER}/{simulation_filename}_{test_case}.pickle', 'rb') as f:
        SS = pickle.load(f)
    return SS


def save_simulation_data(test_case, simulation_filename, data):
    """
    This function saves simulation data to a file.

    Args:
    test_case (str): A string specifying which test case to save.
    simulation_filename (str): A string specifying the name of the simulation.
    data (dict): A dictionary containing the results of the simulation.
    """
    filename = f'{SIMULATIONS_FILE_FOLDER}/{simulation_filename}_{test_case}.pickle'
    filename_smaller = f'{SIMULATIONS_FILE_FOLDER}/{simulation_filename}_{test_case}_smaller.pickle'

    if os.path.exists(filename):
        with open(filename, 'rb') as f:
            old_data = pickle.load(f)
        old_data.update(data)
        data = old_data

    # Save the full version
    with open(filename, 'wb') as f:
        pickle.dump(data, f)

    # Create a smaller version by removing certain keys
    smaller_data = data.copy()
    if 'likelihoods' in smaller_data:
        del smaller_data['likelihoods']
    if 'simulated_chromosomes' in smaller_data:
        del smaller_data['simulated_chromosomes']

    # Save the smaller version
    with open(filename_smaller, 'wb') as f:
        pickle.dump(smaller_data, f)