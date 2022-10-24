## Source Files

**chrom_effects_viz.R** is a R script that takes the m2effects file from simulation outputs and plots where m2 mutations rise in the genome, their effect direction (positive or negative), and their effect size for each phenotype. It writes out plots for the entire genome as well as for each chromosome where m2 mutations arise.

For example, here is a whole genome figure:
![""](/figures/1529274835622_mut_effect_allchrom.png)

And a chromosome figure:
![""](/figures/1529274835622_mut_effect_chrom3.png)

**effects_proc_script.R** is a R script that takes the m2effects file from the simulation outputs and plots the distribution of effect sizes for phenotype0 vs phenotype1. The plot also includes whether the mutation is fixed or not, the correlation for mutations that fixed and those that didn't, and the variance-covariance matrix provided to the simulation.

For example:

![""](/figures/1525515978592_dist_effect_sizes.png)

**phen_env_proc_script.R** is a R script that takes the PhenEnv files from the simulation outputs to understand whether the individuals are actually adapting to the environment by plotting the correlation between the phenotypes and environments that the phenotypes are adapting to for individuals over the simulation generations. If individuals are adapting, the correlation is expected to be high. Individuals are not set to adapt until 4Ne generations so that they can generate genetic variation.

Example for two simulations:

![""](/figures/20180618_phen_env_corr.png)

**interpolating_adaptree_env.R** is an edited R script that takes the original Adaptree data and interpolates the environments into the csv files that can be used as maps in SLiM.

**rep_slim.sh** is a bash script that runs simulations for all combinations of parameter values specified in the script. This script can be modified to include more parameters by adding more for loops.

**slim_sim_count_check.sh** is a bash script that calculates the total number of simulations that will be run given all combinations of parameters values supplied in the script. This helps to determine how many processes will be run on the cluster before running rep_slim.sh.


