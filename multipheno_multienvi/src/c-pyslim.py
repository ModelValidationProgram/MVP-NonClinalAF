import pyslim 
import msprime
import argparse
import numpy as np
import random
import time
import re
import ast
import sys
import os

import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
# parser.add_argument("-s", "--seed", required = True)

args = vars(parser.parse_args())

seed = 892657863
r = 1e-05
N = 1000
mu = 1e-07
output = "multipheno_multienvi/output_multisim/"


# Load the .trees file
T2 = pyslim.load(output+"892657863.trees")

# Calculate tree heights, giving uncoalesced sites the maximum time
def tree_heights(ts):
    heights = np.zeros(ts.num_trees + 1)
    for tree in ts.trees():
        if tree.num_roots > 1: # not fully coalesced
            heights[tree.index] = ts.slim_generation
        else:
            children = tree.children(tree.root)
            real_root = tree.root if len(children) > 1 else children[0]
            heights[tree.index] = tree.time(real_root)
    heights[-1] = heights[-2] # repeat the last entry for plotting
    return heights
 
# Recapitate!
start = time.time()
recapT2 = T2.recapitate(recombination_rate = r, Ne=N, random_seed=seed)
end = time.time()
print("Time it took to to recapitate:", end - start)


# Mutation! 
mutatedT2 = pyslim.SlimTreeSequence(msprime.mutate(recapT2, rate=mu, random_seed=seed, keep=True))

print(f"The tree sequence now has {mutatedT2.num_trees} trees,"
      f" and {mutatedT2.num_mutations} mutations.")
      
# mutated.dump("results/2388682558411_recap_mutate.trees")

# Plot tree hieghts after recapitation
# breakpoints = list(recap.breakpoints())
# heights = tree_heights(recap)
# plt.step(breakpoints, heights, where='post')

# Plot tree heights before recapitation
# breakpoints1 = list(ts.breakpoints())
# heights1 = tree_heights(ts)
# plt.step(breakpoints1, heights1, where='post')
# plt.show()
    

with open(output + str(seed) + "_PlusNeuts.vcf", "w") as vcf_file:
	mutatedT2.write_vcf(vcf_file)