import numpy as np
import sys

# Extracted code to generate the genome doubling Markov matrix G
max_CN = int(sys.argv[1])  # Assuming max_CN is passed as a command-line argument

G = matrix(QQ, max_CN+1, max_CN+1, 0)  # Create a zero matrix of the appropriate size
for i in range(round(max_CN/2)):
    G[i, 2*i] = 1

# Save G to an .sobj file
save(G, 'GD.sobj')

