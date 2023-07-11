import logging
import numpy as np


def get_branch_lengths(trees_and_timings, max_epoch):
    tree, labelled_tree, count, epochs_created, parents = trees_and_timings
    branch_lengths = calculate_child_parent_diff(
        epochs_created = epochs_created, 
        parents = parents,
        max_epoch = max_epoch
    )
    CNs = extract_copy_numbers(tree)
    stacked_branch_lengths,unique_CNs = stack_same_CN_branch_lengths(CNs, branch_lengths)
   
    logging.debug("trees_and_timings")
    logging.debug(str(trees_and_timings))
    logging.debug("max_epoch")
    logging.debug(max_epoch)
    logging.debug("branch_lengths")
    logging.debug(branch_lengths)
    logging.debug("stacked_branch_lengths")
    logging.debug(stacked_branch_lengths)

    return CNs, unique_CNs, branch_lengths, stacked_branch_lengths


def get_path_code(code_list):
    output = ""
    count = 0

    for code in code_list:
        if code == "A":
            count += 1
        elif code == "GD":
            output += str(count)
            count = 0
            output += "G"

    output += str(count)
    return output

# Configure logging settings (you only need to do this once in your script or module)
# this would be a good idea to use throughout the script

def timing_struct_to_all_structures(trees_and_timings, pre, mid, post, max_epoch):
    all_structures = {}
    
    for chrom in trees_and_timings:
        all_structures[chrom] = timing_structs_to_all_structs_per_chrom(trees_and_timings[chrom], pre, mid, post, max_epoch)

    

    return all_structures


def timing_structs_to_all_structs_per_chrom(trees_and_timings, pre, mid, post, max_epoch):
    logging.debug("trees_and_timings")
    logging.debug(trees_and_timings)
    all_structures = [] 
    for index, these_tts in enumerate(trees_and_timings):  # these tts are a 2d array
        if None in these_tts[3]:
            continue
            #BP_likelihoods = -1
        else:
            # trace back to here, asdfasdf
            CNs, unique_CNs, branch_lengths, stacked_branch_lengths = get_branch_lengths(trees_and_timings=these_tts, max_epoch=max_epoch)

            logging.debug("CNs, unique_CNs, branch_lengths, stacked_branch_lengths")
            logging.debug(f"{CNs}, {unique_CNs}, {branch_lengths}, {stacked_branch_lengths}")
            logging.debug("starts and ends")

            path = create_path(pre, mid, post)

            logging.debug(path)

            starts = these_tts[3] #+1
            ends = these_tts[3] + branch_lengths #+1

            logging.debug("starts")
            logging.debug(starts)
            logging.debug("ends")
            logging.debug(ends)

            paths = calculate_BP_paths(branch_lengths, starts, ends, path)

        tree, labelled_tree, count, epochs_created, parents = these_tts

        all_structures += [{
            "pre": pre,
            "mid": mid,
            "post": post,
            "path": path,
            "tree": tree,
            "parents": parents,
            "labelled_tree": labelled_tree,
            "count": count,
            "epochs_created": epochs_created,
            "CNs": CNs,
            "branch_lengths": branch_lengths,
            "unique_CNs": unique_CNs,
            "stacked_branch_lengths": stacked_branch_lengths,
            "starts":starts,
            "ends":ends,
            "paths":paths
        }]

    return all_structures



def timing_struct_to_BP_likelihood_per_chrom(data, structures, p_up, p_down):
    #logging.getLogger().setLevel(logging.DEBUG)

    all_BP_likelihoods = []

    for structure in structures: 
        assert(None not in structure["epochs_created"])
        logging.debug("structures")
        logging.debug(structures)

        likelihoods = calculate_likelihoods_from_paths(paths=structure["paths"], CNs=structure["CNs"], data=data, p_up=p_up, p_down=p_down)

        logging.debug("paths")
        logging.debug(structure["paths"])

        logging.debug("likelihoods")
        logging.debug(likelihoods)

        ll = np.log(likelihoods)
        
        logging.debug("log likelihoods")
        logging.debug(ll)
        structure["BP_individual_log_likelihoods"] = ll

        BP_likelihood = np.sum(ll[:, 1:], axis=1)

        logging.debug("sum them to get BP likelihoods")
        logging.debug(BP_likelihood)

        structure["BP_loglikelihood"] = BP_likelihood

    #logging.getLogger().setLevel(logging.INFO)




def create_path(pre, mid, post):
    path = []
    if pre > 0:
        path += ["A"] * pre
    if mid > -1:
        path += ["GD"]
    if mid > 0:
        path += ["A"] * mid
    if post > -1:
        path += ["GD"]
    if post > 0:
        path += ["A"] * post

    return path


def calculate_BP_paths(branch_lengths, starts, ends, path):
    paths = np.zeros(ends.shape, dtype=object, order="C")

    for row in range(branch_lengths.shape[0]):
        for col in range(branch_lengths.shape[1]):

            these_paths = path[starts[row][col] : ends[row][col]]  # MODIFIED THIS, BUT NOT ENTIRELY SURE, CHECK, ERROR
            path_code = get_path_code(these_paths)

            paths[row][col] = path_code

    # this needs to be updated for paths with length longer than 1. 
    # suppose a path was like GG1. The branching process can take this and calculate a probability, however when the SNV likelihood doe sthe same thing it gives an incompatible probability. 
    # the SNV likelihood needs to know how long SNV can accumulate on a branch. for SNVs to accurately take iunto consideration of how lon g they can accumulate, the BP lieklihood needs to recognise this. 

    # therefore GG1 really needs to look like p(GG1) = p(GD, GD) * p(D)^3 * 4 * p(L)
    # this is not straightforward to program. How do we do it?
    # furthermroe, how do we recognise it?
    #

    # first we recoginise this phenomenum as occuring to paths greater than length 1. all other paths naturally work.
    # P(GG1 1 to 2 true) = P(GG1 1 to 1) * UA/(1-U-D)
    # this works because there si always one chromosome where nothing happens to it. in the 1 to 1 caseA

    # does this also work fro paths of length 1?
    # yes if it is a non genome doubling one. 
    # actually it always works unless the last epoch is a gd one. then it is just the usual. /calculate_BP_paths

    return paths


def check_last_character(code):
    """
    Checks if the last character of the given code is 'G' or '0' versus a non-zero integer.

    Parameters:
        code (str): Code represented as a string.

    Returns:
        bool: True if the last character of the code is 'G' or '0', False if it's a non-zero integer.
    """

    while code and code[-1] == '0':
        code = code[:-1]
        if len(code) == 0: 
            return True

    if code and code[-1] in {'G', '0'}:
        return True
    else:
        return False

def shorten_by_one(path):
    """
    Shortens the given path by one. If the last character is '0', it is removed along with the preceding 'G'. 
    If the path is '0', a ValueError is raised.
    """
    if path == '0':
        raise ValueError("Cannot shorten '0'")
    
    if path.endswith('G0'):
        return path[:-2]
    elif path[-1].isdigit() and int(path[-1]) > 0:
        return path[:-1] + str(int(path[-1]) - 1)
    else:
        raise ValueError(f"Invalid path: {path}")


def calculate_likelihoods_from_paths(paths, CNs, data, p_up, p_down):
    likelihoods = np.zeros(paths.shape, dtype=float, order="C")

    assert isinstance(p_up, int) and 0 <= p_up <= 100, f"p_up should be an integer from 0 to 100: {p_up}"
    assert isinstance(p_down, int) and 0 <= p_down <= 100, f"p_down should be an integer from 0 to 100: {p_down}"
    p_up = p_up/100.0
    p_down = p_down/100.0

    assert 0 <= p_up <= 1 and 0 <= p_down <= 1, f"p_up: {p_up}, p_down: {p_down}, p_up and p_down should be between 0 and 1"


    for row in range(paths.shape[0]):
        for col in range(paths.shape[1]):
            path = paths[row][col]
            if check_last_character(path) or CNs[col] < 2:
                likelihood = data[path][1][min(2,CNs[col])] # the prob of going from CN 1 to ...
                # the reason for the min(2,*) is that if it is a copy number zero node then find the prob it is zero, if it is 1 then likewise, if 2 or more then it is just the prob of a bifuctation\
            else:
                # for the sake of the SNV likelihood we need the doubling up to be at the last second:
                if len(path) > 1:
                    likelihood = data[shorten_by_one(path)][1][1] * p_up
                    n = data[shorten_by_one(path)].shape[1]  # assuming 'data' is a numpy array

                    for j in range(2, n):
                        likelihood += data[shorten_by_one(path)][1][j] * p_up * (p_down ** (j-1)) * j

                else:
                    assert(len(path) == 1)
                    likelihood = data[path][1][2]



            if math.isnan(likelihood):
                logging.getLogger().setLevel(logging.DEBUG)
                logging.debug("Last element in path: %s", path[-1])
                logging.debug("Value of CNs[col]: %s", CNs[col])
                logging.debug("Value of p_up: %s", p_up)
                logging.debug("Value of p_down: %s", p_down)
                logging.debug("Calculated likelihood: %s", likelihood)
                logging.debug("Exiting the program.")
                logging.debug("Likelihood contains 'nan'. Exiting the program.")
                sys.exit()

            likelihoods[row][col] = likelihood


    if np.isnan(likelihoods).any():
        print("The likelihoods array contains NaN values. Exiting the program.")
        sys.exit()

    return likelihoods



def add_BP_likelihoods(all_structures, p_up, p_down):

    file = PRECOMPUTED_FILE_FOLDER + "/subbed_mat_u"+str(int(p_up))+"_d"+str(int(p_down))+"_p8_v4.precomputed_paths.pickle"

    data = pkl.load(open(file,'rb'))


    for chrom in all_structures.keys():
        timing_struct_to_BP_likelihood_per_chrom(
            data=data,
            structures=all_structures[chrom],
            p_up=p_up,
            p_down=p_down
            )
    return all_structures