from typing import Dict



def analyze_sorted_results(SS):
    print("Do some analysis")
    for solution in SS["solutions"]:
        total_truth_nodes = 0
        total_solution_nodes = 0
        num_chrom_with_correct_CN = 0
        num_chrom_with_correct_CN_and_epoch_created = 0
        average_distance_from_truth_of_epoch_created = 0
        SS["total_nodes"] = {}

        for chrom in range(23):
            # are trees same by epoch and copy number


            # count the total nodes:
            solution[chrom]["total_nodes"] = count_nodes(solution[chrom]["dict_tree"])
            total_solution_nodes += solution[chrom]["total_nodes"]

            SS["total_nodes"][chrom] = count_nodes(SS["simplified_truth_trees"][chrom])
            total_truth_nodes += SS["total_nodes"][chrom]

            # count the number with the correct CN:
            solution[chrom]["correct_CN"] = analyze_tree_copy_number(SS["simplified_truth_trees"][chrom], solution[chrom]["dict_tree"])
            num_chrom_with_correct_CN += solution[chrom]["correct_CN"]

            # count the number with the correct CN and epoch created:
            solution[chrom]["correct_CN_and_epoch_created"] = analyze_tree_properties(SS["simplified_truth_trees"][chrom], solution[chrom]["dict_tree"])
            num_chrom_with_correct_CN_and_epoch_created += solution[chrom]["correct_CN_and_epoch_created"]

            # find the average epoch distance for tree of correct CN (incorrect CN trees get infinite distance)
            solution[chrom]["error_in_epoch_created_estimate"] = calculate_tree_distance(SS["simplified_truth_trees"][chrom], solution[chrom]["dict_tree"], total_nodes)
            average_distance_from_truth_of_epoch_created += solution[chrom]["error_in_epoch_created_estimate"]


    print_summary(total_nodes, num_chrom_with_correct_CN,
                  num_chrom_with_correct_CN_and_epoch_created, average_distance_from_truth_of_epoch_created)

    for result_index, res in enumerate(sorted_results):
        new_result = compile_new_result(res, SS, pre_est, mid_est, post_est)

        process_new_result(new_result, SS)

        file_name = f'{SIMULATIONS_FILE_FOLDER}/{simulation_filename}_{test_case}.pickle'
        all_results = load_simulation_data(test_case, simulation_filename)
        all_results[str(result_index)] = new_result
        write_simulation_data(test_case, simulation_filename, all_results)

def get_all_trees_ready_for_comparison(SS):
    #print("get trees ready for comparison")

    get_truth_trees_ready_for_comparison(SS)
    #print("truth trees are ready")

    get_estimated_trees_ready_for_comparison(SS)
    #print("estimated trees are ready")


def get_truth_trees_ready_for_comparison(SS):
    #pretty_print("simplified simulated tree")

    if "simplified_truth_trees" not in SS:
        SS["simplified_truth_trees"] = {}

    assert(sorted(list(SS["truth_trees"].keys())) == list(range(23))) # assert it has all the keys available

    for chrom in SS["truth_trees"]:
        #print("Current chromosome: ", chrom)
        SS["simplified_truth_trees"][chrom] = filter_tree(SS["truth_trees"][chrom],keys_to_keep = ["copy_number","epoch_created"])
        #print("Newly added tree to simplified_truth_trees: ", SS["simplified_truth_trees"][chrom])

        SS["simplified_truth_trees"][chrom] = create_epoch_index(SS["simplified_truth_trees"][chrom], key_from="epoch_created", key_to="epoch_index")
        #print("Trees after adding epoch index: ",  SS["simplified_truth_trees"][chrom])

        SS["simplified_truth_trees"][chrom] = order_tree_keys_alphabetically(SS["simplified_truth_trees"][chrom])
        #print("Tree after ordering keys alphabetically: ", SS["simplified_truth_trees"][chrom])

        SS["simplified_truth_trees"][chrom] = sort_tree_by_copynumber(SS["simplified_truth_trees"][chrom])
        #print("Tree after ordering keys by copynumber: ", SS["simplified_truth_trees"][chrom])
        #print("\n"*5)


def get_estimated_trees_ready_for_comparison(SS):
    for solution in SS["solutions"]:
        for chrom in range(23):
            assert(chrom in solution)
            solution[chrom]["dict_tree"] = convert_CN_tree_and_epoch_list_to_dict_tree(solution[chrom]["tree"], solution[chrom]["epochs_created"])

            solution[chrom]["dict_tree"] = create_epoch_index(solution[chrom]["dict_tree"], key_from="epoch_created", key_to="epoch_index")

            solution[chrom]["dict_tree"] = order_tree_keys_alphabetically(solution[chrom]["dict_tree"])

            solution[chrom]["dict_tree"] = sort_tree_by_copynumber(solution[chrom]["dict_tree"])

            solution[chrom]["dict_tree"]['epoch_index'] = 0


def sort_simulation_results_by_likelihood(solutions):
    """
    This function sorts simulation results based on 'best_loglik'.
    
    Args:
    results (list): A list of dictionaries containing simulation results.

    Returns:
    list: A sorted list of dictionaries.
    """
    return sorted(solutions, key=lambda x: x['est_neg_loglik'], reverse=True)






def add_ev_strings_and_counts_to_dicts(SS):
    truth_str = get_ev_string(SS['pre'], SS['mid'], SS['post'])
    true_GD_count = count_genome_doublings(truth_str)
    for d in SS["solutions"]:
        d['ev_str'] = truth_str
        d['est_pre'] = d[0]['pre']
        d['est_mid'] = d[0]['mid']
        d['est_post'] = d[0]['post']
        d['pre'] = SS['pre']
        d['mid'] = SS['mid']
        d['post'] = SS['post']
        d['p_up'] = SS['p_up']
        d['p_down'] = SS['p_down']
        d['plambda'] = SS['rate']

        d['est_ev_str'] = get_ev_string(d['est_pre'], d['est_mid'], d['est_post'])
        d['genome_doublings'] = true_GD_count
        d['est_genome_doublings'] = count_genome_doublings(d['est_ev_str'])


def add_aic_to_dicts(SS):
    for d in SS["solutions"]:
        num_parameters = d['est_genome_doublings'] + 3
        neg_log_likelihood = d['est_neg_loglik']
        aic = 2 * num_parameters + 2 * neg_log_likelihood
        d['AIC'] = aic


def sort_dicts_by_worst_aic(SS):
    SS["solutions"] = sorted(SS["solutions"], key=lambda x: -x['AIC'])

def count_genome_doublings(ev_string):
    return ev_string.count("G")

def process_further(SS):
    add_ev_strings_and_counts_to_dicts(SS)
    add_aic_to_dicts(SS)
    sort_dicts_by_worst_aic(SS)

