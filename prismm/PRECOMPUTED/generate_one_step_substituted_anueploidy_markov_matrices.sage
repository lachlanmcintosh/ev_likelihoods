# generate_substituted_matrices.sage
import sys
import os
import pickle
from sage.all import *

args = sys.argv
p_up = args[1]
p_down = args[2]
path_description = args[3]

# Load the anueploidy Markov matrix
anueploidy_matrix = load("MATRICES/matrix_"+path_description+".sobj")

# Substitute the up and down values
subbed_matrix = anueploidy_matrix.subs(u=int(p_up)/100, d=int(p_down)/100)

# Convert to double precision and save as pickle
subbed_matrix = subbed_matrix.apply_map(RR)
subbed_matrix = subbed_matrix.numpy(dtype='double')

base_filename = "MATRICES/subbed_mat_u"+p_up+"_d"+p_down+"_"+path_description+".pickle"
with open(base_filename, 'wb') as m_output:
    pickle.dump(subbed_matrix, m_output)

