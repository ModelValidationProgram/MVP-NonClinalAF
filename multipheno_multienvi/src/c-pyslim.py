import msprime, pyslim, argparse

import numpy as np

import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--seed", required = True)

args = vars(parser.parse_args())

seed = 892657863
r = 1e-05
N = 1000
mu = 1e-07


# Load the .trees file
T2 = pyslim.load(seed + "_Invers.trees")

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
recapT2 = T2.recapitate(recombination_rate = r, Ne=N, random_seed=seed)


# Mutation! 
mutatedT2 = pyslim.SlimTreeSequence(msprime.mutate(recapT2, rate=mu, random_seed=seed, keep=True))

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