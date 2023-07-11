import json
import sys
import pandas as pd
from typing import Dict


def pretty_print_tree(tree):
    """
    Pretty-prints a nested dictionary representing a tree.

    Args:
        tree (dict): The tree dictionary to be pretty-printed.
    """
    print(json.dumps(tree, indent=4, sort_keys=True))



def print_summary(total_nodes, num_chrom_with_correct_CN, num_chrom_with_correct_CN_and_epoch_created, average_distance_from_truth_of_epoch_created):
    print(f"Total nodes: {total_nodes}")
    print(f"Number of chromosomes with correct copy numbers: {num_chrom_with_correct_CN}")
    print(f"Number of chromosomes with correct copy numbers and epochs created: {num_chrom_with_correct_CN_and_epoch_created}")
    print(f"Average distance from truth of epoch created: {average_distance_from_truth_of_epoch_created}")


def print_tree(tree):
    """
    Takes a dictionary tree and prints a neatly formatted list representation.
    
    :param tree: Dictionary tree to be printed.
    """
    converted_tree = convert_dict_tree_to_list(tree)
    
    def print_tree_recursively(tree, indent=0):
        """
        Helper function to recursively print a tree.
        
        :param tree: List representation of the tree.
        :param indent: Current indentation level (default: 0).
        """
        print('\t' * indent, tree[0])
        for subtree in tree[1:]:
            if subtree is not None:
                print_tree_recursively(subtree, indent + 1)

    print_tree_recursively(converted_tree)


def print_solution_results(SS):
    """
    Prints the solution results for intuition. For each chromosome in the solution,
    this function prints the dictionary representation of the solution tree and the corresponding
    truth tree from the simplified truth trees.

    :param SS: A dictionary containing the solutions and the simplified truth trees.
    """
    print("Print results for intuition:")
    for idx, solution in enumerate(SS["solutions"]):
        print(f"\nSolution {idx}:")

        # Print non-integer keys:
        for key, value in solution.items():
            if not isinstance(key, int):
                print(f"{key}: {value}")

        for chrom in range(23):
            print("\n\tSolution Tree for Chromosome", chrom)
            print(solution[chrom]["dict_tree"])
            print_tree(solution[chrom]["dict_tree"])
            print("\n\tTruth Tree for Chromosome", chrom)
            print(SS["simplified_truth_trees"][chrom])
            print_tree(SS["simplified_truth_trees"][chrom])



def print_similarity_results(SS: Dict):
    """
    Print the similarity results computed in each solution along with comparison results of relative timings.

    :param SS: A dictionary containing the solutions, the computed similarity results, and comparison results.
    """
    print("Print similarity results:")

    # Set the display options.
    pd.set_option('display.max_columns', None) 
    pd.set_option('display.max_rows', None) 
    pd.set_option('display.max_colwidth', None)
    pd.set_option('display.width', 200)  # default is 80


    for idx, solution in enumerate(SS["solutions"]):
        print(f"\nSolution {idx}:")

        # Initialize an empty DataFrame for each solution
        df = pd.DataFrame(columns=['Copy Number Distance Function', 'Epoch Created Distance Function', 'Normalising Constant', 'Similarity Score'])

        # Print the keys related to similarity results
        for key, value in solution.items():
            if isinstance(key, str):
                if any(distance_function in key for distance_function in ["ZeroDistance", "AbsoluteDifference", "SquaredDistance", "InfiniteDistance"]):
                    distance_function_1, distance_function_2, constant = key.split('_')

                    # Append a new row to the DataFrame
                    new_row = pd.DataFrame([{
                        'Copy Number Distance Function': distance_function_1,
                        'Epoch Created Distance Function': distance_function_2,
                        'Normalising Constant': constant,
                        'Similarity Score': value
                    }])
                    df = pd.concat([df, new_row], ignore_index=True)

        solution["distance_dataframe"] = df

        # Print the DataFrame
        df = df.sort_values(by=['Normalising Constant', 'Copy Number Distance Function'])
        print(df)

        solution["correct_path"] = SS['pre'] == solution[0]['pre'] and SS['mid'] == solution[0]['mid'] and SS['post'] == solution[0]['post']
         
        solution_data = {
            'best_neg_loglik': [solution['est_neg_loglik']],
            'best_p_up': [solution['est_p_up']],
            'best_p_down': [solution['est_p_down']],
            'best_plambda': [solution['est_plambda']],
            'pre': [solution[0]['pre']],
            'mid': [solution[0]['mid']],
            'post': [solution[0]['post']],
            'total_time': [solution['execution_times']['total_time']],
            'time_get_all_trees_and_timings': [solution['execution_times']['get_all_trees_and_timings']],
            'time_timing_struct_to_all_structures': [solution['execution_times']['timing_struct_to_all_structures']],
            'time_find_BP_and_SNV_loglik': [solution['execution_times']['timing_struct_to_all_structures']],
            'counts_<': [solution['counts_<']],
            'counts_<=': [solution['counts_<=']],
            'correct_path': solution["correct_path"]
        }

        df = pd.DataFrame(solution_data)
        print(df)
        solution["metadata_dataframe"] = df



def get_size(obj):
    size = sys.getsizeof(obj)
    if isinstance(obj, dict):
        size += sum(get_size(value) for value in obj.values())
    elif isinstance(obj, list):
        size += sum(get_size(item) for item in obj)
    return size

def print_dicts_memory(ss, indent=0, level=1, max_level=4):
    if level > max_level:
        return

    for key, value in ss.items():
        size = get_size(value)
        if isinstance(value, dict):
            print(' ' * indent + str(key) + ': dict of total size ' + str(size) + ' bytes')
            print_dicts_memory(value, indent + 4, level + 1, max_level)
        elif isinstance(value, list):
            print(' ' * indent + str(key) + ': list of total size ' + str(size) + ' bytes')
        else:
            print(' ' * indent + str(key) + ': ' + str(size) + ' bytes')

