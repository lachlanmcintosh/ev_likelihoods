from collections import defaultdict
import numpy as np
import pandas as pd


###### STEP 3; calculate the SNV multiplicities of each chromosome

def simulated_chromosomes_to_SNV_counts(simulated_chromosomes):
    """
    Count the copy number of each SNV in the simulated chromosomes.

    :param simulated_chromosomes: A dictionary containing the simulated chromosomes.
    :return: A dictionary with the count of copy numbers for each SNV in the simulated chromosomes.
    """

    SNV_copy_counter = {}

    for chrom_type in simulated_chromosomes:
        if chrom_type not in SNV_copy_counter:
            SNV_copy_counter[chrom_type] = {}

        for chrom in simulated_chromosomes[chrom_type]:
            if chrom["dead"]:
                continue

            for SNV in chrom["SNVs"]:
                UI = SNV["unique_identifier"]

                if UI not in SNV_copy_counter[chrom_type]:
                    SNV_copy_counter[chrom_type][UI] = 1
                else:
                    SNV_copy_counter[chrom_type][UI] += 1

    return SNV_copy_counter



def SNV_counts_to_SNV_multiplicities(SNV_copy_counter):
    """
    Convert the copy number counts of each SNV into SNV multiplicities.

    :param SNV_copy_counter: A dictionary with the count of copy numbers for each SNV.
    :return: A dictionary with the SNV multiplicities.
    """

    multiplicities = {}

    for chrom_number in SNV_copy_counter:
        if chrom_number not in multiplicities:
            multiplicities[chrom_number] = {}

        for SNV in SNV_copy_counter[chrom_number]:
            CN = SNV_copy_counter[chrom_number][SNV]

            if CN not in multiplicities[chrom_number]:
                multiplicities[chrom_number][CN] = 1
            else:
                multiplicities[chrom_number][CN] += 1

    return multiplicities


def count_SNV_multiplicities(simulated_chromosomes):
    return SNV_counts_to_SNV_multiplicities(simulated_chromosomes_to_SNV_counts(simulated_chromosomes))



