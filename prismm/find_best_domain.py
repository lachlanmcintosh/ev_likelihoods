from datetime import datetime, timezone
import os
import glob
import pickle
import pandas as pd

def load_results(filename):
    with open(filename, 'rb') as f:
        data = pickle.load(f)
    return data

def get_values(data):
    proportion_estimated = data.get('proportion_estimated', None)
    best_accuracy = data.get('best_accuracy', None)
    best_estimate_accuracy = data.get('best_estimate_accuracy', None)
    arbitrary_estimate_accuracy = data.get('arbitrary_estimate_accuracy', None)
    return proportion_estimated, best_accuracy, best_estimate_accuracy, arbitrary_estimate_accuracy


def get_sorted_filenames(directory, prefix):
    files = glob.glob(os.path.join(directory, f'{prefix}_lam*_alpha*_summary.pkl'))
    special_file = os.path.join(directory, f'{prefix}_summary.pkl')
    if os.path.exists(special_file):
        files.append(special_file)

    data = []
    for file in files:
        # skip files where 'lam' is 0.5
        if 'lam0.5' in file:
            continue
        results = load_results(file)
        proportion_estimated, best_accuracy, best_estimate_accuracy, arbitrary_estimate_accuracy = get_values(results)
        if proportion_estimated is not None and best_accuracy is not None:
            diff = proportion_estimated - best_accuracy
            is_average_estimate = file == special_file
            data.append((file, proportion_estimated, best_accuracy, best_estimate_accuracy, arbitrary_estimate_accuracy, diff, is_average_estimate))
    data.sort(key=lambda x: x[5], reverse=True)
    return data


def print_sorted_files(directory, prefix):
    data = get_sorted_filenames(directory, prefix)
    for file, proportion_estimated, best_accuracy, best_estimate_accuracy, arbitrary_estimate_accuracy, diff, is_average_estimate in data:
        print(f"Filename: {file}\nProportion_estimated: {proportion_estimated}\nBest accuracy: {best_accuracy}\nBest estimate accuracy: {best_estimate_accuracy}\nArbitrary estimate accuracy: {arbitrary_estimate_accuracy}\nDifference: {diff}\nAverage Estimate: {is_average_estimate}\n")

def create_dataframe(input_directory, output_directory, prefix):
    data = get_sorted_filenames(input_directory, prefix)
    df = pd.DataFrame(data, columns=["filename", "proportion_estimated", "best_accuracy", "best_estimate_accuracy", "arbitrary_estimate_accuracy", "difference", "average_estimate"])
    df.to_csv(os.path.join(output_directory, 'summary.csv'), index=False)

import argparse

def main(input_directory, output_directory, prefix):
    print_sorted_files(input_directory, prefix)
    create_dataframe(input_directory, output_directory, prefix)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_directory", default='SIMULATIONS')
    parser.add_argument("--output_directory", default="~/ev_likelihood/")
    parser.add_argument("--prefix", default='simulation')
    args = parser.parse_args()
    main(args.input_directory, args.output_directory, args.prefix)

