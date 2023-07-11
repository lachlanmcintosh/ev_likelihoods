# constants.py

BASE_FILE_FOLDER = "/vast/scratch/users/lmcintosh/CN_SV_SNV_evolutionary_likelihood/"
PRECOMPUTED_FILE_FOLDER = BASE_FILE_FOLDER + "PRECOMPUTED/MATRICES/"
SIMULATIONS_FILE_FOLDER = BASE_FILE_FOLDER + "SIMULATIONS"


# Percent of total (Female) genome
LENGTHS = {
    0: 8.18,
    1: 8.04,
    2: 8.60,
    3: 6.33,
    4: 5.98,
    5: 5.65,
    6: 5.25,
    7: 4.84,
    8: 4.64,
    9: 4.48,
    10: 4.45,
    11: 4.38,
    12: 3.78,
    13: 3.52,
    14: 3.32,
    15: 2.94,
    16: 2.61,
    17: 2.52,
    18: 2.11,
    19: 2.07,
    20: 1.55,
    21: 1.64,
    22: 5.13
}

# adjust the rate of how frequently SNVs occur per genome, so that longer chromosomes have a proportionally larger chance of a new SNV
for chrom_type in LENGTHS:
    LENGTHS[chrom_type] *= 1 / 100


