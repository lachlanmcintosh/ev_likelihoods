import numpy as np
import scipy.special
import sys
import logging


def get_poisson_loglikelihood(counts, stacked_branch_lengths, plambda, chrom, lengths, unique_CNs):

    stacked_branch_lengths_float = stacked_branch_lengths.astype(float)

    A = np.log(stacked_branch_lengths_float * plambda * lengths[chrom]) * counts
    A = np.where((np.isnan(A)) & (stacked_branch_lengths_float == 0), -np.inf, A)

    B = -np.tile([scipy.special.gammaln(x + 1) for x in counts], (stacked_branch_lengths.shape[0], 1))

    C = -stacked_branch_lengths_float * plambda * lengths[chrom]

    summed = A + B + C
    total = np.sum(summed, axis=1)

    summed_numeric = summed.astype(float)

    if np.isnan(summed_numeric.flatten()).any():
        logging.warning("#####################\n"*10)
        logging.warning("counts: %s", counts)
        logging.warning("plambda: %s", plambda)
        logging.warning("chrom: %s", chrom)
        logging.warning("lengths: %s", lengths)

        for i in range(stacked_branch_lengths.shape[0]):
            logging.warning("Row %s:", i+1)
            logging.warning("stacked_branch_lengths: %s", stacked_branch_lengths[i])
            logging.warning("A: %s", A[i])
            logging.warning("B: %s", B[i])
            logging.warning("C: %s", C[i])
            logging.warning("summed: %s", summed[i])

        logging.warning("total: %s", total)
        sys.exit()

    # FIX 1, make sure that were not testing if 0 SNVs occured in 0 time:
    # Find the indices where branch_lengths is zero
    zero_indices_branch_lengths = np.argwhere(stacked_branch_lengths == 0)

    # Make sure that the position in summed_numeric corresponding to a zero in branch_lengths and counts is also zero
    for idx in zero_indices_branch_lengths:
        if counts[idx[1]] == 0:
            summed_numeric[idx[0], idx[1]] = 0  # set the specific index to 0


    # FIX 2, make sure that we are not testing to see SNVs on CN 0 chromosomes:
    # If CNs is zero in any column, set every entry in summed_numeric in that column to zero
    for j, cn in enumerate(unique_CNs):
        if cn == 0:
            summed_numeric[:, j] = 0  # set the entire column to 0


    return summed_numeric


def test_get_poisson_loglikelihood():
    counts = np.array([2, 3, 5, 7, 11])
    stacked_branch_lengths = np.array(
        [
            [0.5, 0.6, 0.7, 0.8, 0.9],
            [0.1, 0.2, 0.3, 0.4, 0.5],
        ]
    )
    plambda = 0.2
    chrom = "chrom1"
    to_delete = [1, 3]
    lengths = {"chrom1": 1000}
    expected_total = [np.array([-351.94186892 -125.86288737]), np.array([-477.80475629])]
    total = get_poisson_loglikelihood(counts, stacked_branch_lengths, plambda, chrom, to_delete, lengths)
    assert np.allclose(total, expected_total, rtol=1e-9, atol=1e-9), f"get_poisson_loglikelihood test failed, expected: {expected_total}, but got: {total}"

#test_get_poisson_loglikelihood()

def get_all_poisson_loglikelihoods_per_chr(chrom_structures, plambda, observed_SNV_multiplicities, chrom):
    SNV_likelihoods = []

    for i in range(len(chrom_structures)):
        # get the counts of SNVs per unique CN
        counts = []
        unique_CNs = chrom_structures[i]["unique_CNs"]
        for CN in unique_CNs:
            if int(CN) not in observed_SNV_multiplicities:
                counts += [0]
            else:
                counts += [observed_SNV_multiplicities[CN]]

        assert(unique_CNs == sorted(unique_CNs,reverse=True))

        this_SNV_likelihoods = get_poisson_loglikelihood(
            counts=counts,
            stacked_branch_lengths=chrom_structures[i]["stacked_branch_lengths"],
            plambda=plambda,
            chrom=chrom,
            lengths=LENGTHS,
            unique_CNs=unique_CNs
        )
        
        this_SNV_likelihood = np.sum(this_SNV_likelihoods, axis=1)
        chrom_structures[i]["SNV_log_likelihoods"] = this_SNV_likelihoods #np.transpose(this_SNV_likelihoods)
        #sys.exit()

        chrom_structures[i]["counts"] = counts
        chrom_structures[i]["SNV_log_likelihood"] = this_SNV_likelihood
        chrom_structures[i]["SNV_BP_log_likelihood"] = chrom_structures[i]["SNV_log_likelihood"] + chrom_structures[i]["BP_loglikelihood"]
        SNV_likelihoods += [chrom_structures[i]["SNV_BP_log_likelihood"]]
        chrom_structures[i]["observed_SNV_multiplicities"]  = observed_SNV_multiplicities
    return SNV_likelihoods


def test_get_all_poisson_loglikelihoods_per_chr():
    timings = [
        (None, None, None, [0, 1, 2], None),
        (None, None, None, [0, 2, 4], None)
    ]
    plambda = 0.2
    BP_likelihoods = [np.array([-1, -2, -3]), np.array([-4, -5, -6])]
    observed_SNV_multiplicities = {1: 3, 2: 5}
    chrom = "chrom1"

    expected_SNV_likelihoods = [
        np.array([-17.39516599, -19.73780159, -22.08043718]),
        np.array([-21.13183242, -23.47446761, -25.8171028])
    ]

    SNV_likelihoods = get_all_poisson_loglikelihoods_per_chr(timings, plambda, BP_likelihoods, observed_SNV_multiplicities, chrom)

    for i in range(len(SNV_likelihoods)):
        assert np.allclose(SNV_likelihoods[i], expected_SNV_likelihoods[i], rtol=1e-9, atol=1e-9)



# test_get_all_poisson_loglikelihoods_per_chr()


def find_best_SNV_likelihood(plambda, these_structures, observed_SNV_multiplicities):
    SNV_likelihoods = {}
    best = {}
    chroms = these_structures.keys()

    for chrom in chroms:
        # ERROR SNV_likelihoods[chrom] is saved in these_structures so doesn't need to get passed back like this
        SNV_likelihoods[chrom] = get_all_poisson_loglikelihoods_per_chr(
            chrom_structures = these_structures[chrom],
            plambda = plambda,
            observed_SNV_multiplicities = observed_SNV_multiplicities[chrom],
            chrom=chrom
        )
        
        best[chrom] = (-np.Inf, 0, 0)

        for tree_index in range(len(SNV_likelihoods[chrom])):
            max_SNV_loglik = max(SNV_likelihoods[chrom][tree_index])
            best_timings_row_index = np.argmax(SNV_likelihoods[chrom][tree_index])

            if max_SNV_loglik > best[chrom][0]:
                best[chrom] = (max_SNV_loglik, tree_index, best_timings_row_index)

    total = sum([best[chrom][0] for chrom in chroms])
    return total, best


def test_find_best_SNV_likelihood():
    plambda = 0.2
    timings = {
        "chrom1": [
            (None, None, None, [0, 1, 2], None),
            (None, None, None, [0, 2, 4], None)
        ]
    }
    BP_likelihoods = {
        "chrom1": [np.array([-1, -2, -3]), np.array([-4, -5, -6])]
    }

    expected_total = -17.39516599
    expected_best = {
        "chrom1": (-17.39516599, 0, 0)
    }

    total, best = find_best_SNV_likelihood(plambda, timings, BP_likelihoods,observed_SNV_multiplicities)

    assert np.isclose(total, expected_total, rtol=1e-9, atol=1e-9)

    for chrom in best:
        assert np.isclose(best[chrom][0], expected_best[chrom][0], rtol=1e-9, atol=1e-9)
        assert best[chrom][1:] == expected_best[chrom][1:]