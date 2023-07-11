import copy
import logging
import numpy as np
from scipy import optimize as opt
import sys


def get_best_struct(best_structure_indicies,best_structures):
    best_structure = {}
    for chrom in best_structures:
        best_tree = best_structure_indicies[chrom][1]
        best_structure[chrom] = copy.deepcopy(best_structures[chrom][best_tree])

        best_row = best_structure_indicies[chrom][2]
        for key in ['epochs_created', 'branch_lengths', 'stacked_branch_lengths', 'starts', 'ends', 'paths', 'BP_individual_log_likelihoods', 'BP_loglikelihood', 'SNV_log_likelihood', 'SNV_BP_log_likelihood', 'SNV_log_likelihoods']:
            best_structure[chrom][key] = best_structure[chrom][key][best_row]

    return(best_structure)


def are_all_chromosomes_viable_by_BP_likelihood(all_structures, tree_flexibility):
    """
    Check if all chromosomes have a viable tree.

    A tree is considered viable if its maximum depth is less than or equal to
    the calculated epochs.

    Args:
    all_structures (dict): A dictionary where each key is a chromosome and its
        value is a list of results. Each result is another dictionary that has
        'tree', 'pre', 'mid' and 'post' keys.

    Returns:
    bool: True if all chromosomes have at least one viable tree, False otherwise.
    """
    return_val = True

    for chrom in all_structures:
        tests = [tree_in_bounds(result["tree"], calculate_epochs(result["pre"],result["mid"],result["post"]) ,tree_flexibility) for result in all_structures[chrom]]
        #tests = [max_tree_depth(result["tree"]) <= calculate_max_allowed_BP_depth(result["pre"], result["mid"], result["post"]) for result in all_structures[chrom]]
        if not any(tests):
            return_val = False

    if not return_val:
        logging.info("print what failed")
        logging.info("all_structures %s", all_structures)
        result = all_structures[0][0]
        max_epoch = calculate_epochs(result["pre"],result["mid"],result["post"])
        logging.info("max_epochs %s", max_epoch)
        for chrom in all_structures:
            logging.info("chrom %s", chrom)
            for result in all_structures[chrom]: 
                logging.info("print the trees: %s", result["tree"] )
                logging.info("prihnt the tree depths %s", max_tree_depth(result["tree"]))
            tests = [tree_in_bounds(result["tree"], calculate_epochs(result["pre"],result["mid"],result["post"]) ,tree_flexibility) for result in all_structures[chrom]]
            logging.info("tests %s", tests)
            logging.info("any tests %s", any(tests))
            logging.info("tree_flexibility")
            logging.info(f"DEPTH >= {max_epoch + 3 - tree_flexibility} and depth <= {max_epoch + 3}")


    return return_val



def find_BP_and_SNV_loglik(plambda_start, p_up_start, p_down_start, p_window, plambda_window, all_structures, observed_SNV_multiplicities, total_SNVs, tree_flexibility):
    best_neg_loglik = float("inf")
    best_p_up = 0
    best_p_down = 0
    best_plambda = 0
    best_structures = None
    for p_up in range(max(0,p_up_start-p_window), min(100,p_up_start+p_window+1)): #[p_up_start]:
        for p_down in range(max(0,p_down_start-p_window), min(100,p_down_start+p_window+1)): 
            if p_up+p_down>100:
                continue
            # get the branching process loglik here, 
            # then can pass it in to optimise over the SNVs, otherwise it gets recomputed so many times
            add_BP_likelihoods(all_structures, p_up, p_down)
            temp_structures = copy.deepcopy(all_structures)
            
            def optimize_func(plambda):
                if total_SNVs > 0:
                    assert plambda != 0, f"plambda: {plambda}, observed_SNV_multiplicities: {observed_SNV_multiplicities}"

                total, best = find_best_SNV_likelihood(
                    plambda=plambda,
                    these_structures=all_structures,
                    observed_SNV_multiplicities=observed_SNV_multiplicities
                    )
                return -total

            bounds = (plambda_start * plambda_window, plambda_start / plambda_window)
            x0 = plambda_start
            res = opt.minimize_scalar(optimize_func, bounds=bounds, method='bounded', options={'disp': True})
            #for chrom in all_structures:
            #    for result in all_structures[chrom]:
            #        has_non_inf_values = np.logical_not(np.isinf(result["SNV_BP_log_likelihood"])).any()
            #        if not has_non_inf_values:
            #            print("chrom:", chrom)
            #            print(result)

            if res.fun < best_neg_loglik:
                best_neg_loglik = res.fun
                best_p_up = p_up
                best_p_down = p_down
                best_plambda = res.x
                best_structures = temp_structures #copy.deepcopy(all_structures)
                best_res = res

    # need to work out why this is getting None when it shouldn't
    if best_structures == None:
        # find out why is it none:
        # is it because that for one chromosome each of the tree structures are too deep? check it:
        if not are_all_chromosomes_viable_by_BP_likelihood(all_structures, tree_flexibility):
            print("ERROR, all chroms have viable trees, need to implement code and investigate why the SNVs are throwing this off.")
            sys.exit()
        return None
    total, best_structure_indicies = find_best_SNV_likelihood(
        plambda=best_plambda,
        these_structures=best_structures,
        observed_SNV_multiplicities=observed_SNV_multiplicities
    )
    assert best_neg_loglik == -total
    best_structure = get_best_struct(best_structure_indicies,best_structures)

    return best_neg_loglik, best_p_up, best_p_down, best_plambda, best_structure


def test_find_BP_and_SNV_loglik():
    plambda_start = 0.2
    p_up_start = 50
    p_down_start = 50
    trees_and_timings = {
        "chrom1": [
            (None, None, None, [0, 1, 2], None),
            (None, None, None, [0, 2, 4], None)
        ]
    }
    pre = 0
    mid = 0
    post = 0
    p_window = 0
    plambda_window = 1

    expected_loglik = 17.39516599
    expected_best_p_up = 50
    expected_best_p_down = 50
    expected_best_plambda = 0.2

    loglik, best_p_up, best_p_down, best_plambda, res = find_BP_and_SNV_loglik(plambda_start, p_up_start, p_down_start, trees_and_timings, pre, mid, post, p_window, plambda_window, observed_SNV_multiplicities)

    assert np.isclose(loglik, expected_loglik, rtol=1e-9, atol=1e-9)
    assert best_p_up == expected_best_p_up
    assert best_p_down == expected_best_p_down
    assert np.isclose(best_plambda, expected_best_plambda, rtol=1e-9, atol=1e-9)

