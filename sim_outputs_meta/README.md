Below is a description of metadata for each set of output tables. The first number is the seed for the simulation.

### 2549986039929_fitnessmat_ind.txt
* an n_deme x n_individual table that indicates the fitness of each individual (in columns) in every metapopulation deme (in rows)

### 2549986039929_fitnessmat.txt
* an n_deme x n_deme table that indicates the mean fitness of individuals from the source deme (in columns) in the transplant deme (in rows) CHECK THIS IS CORRECT

### 2549986039929_ind.txt
Individual information output last generation

* `seed` simulation seed
* `indID` individual ID in SliM. Corresponds to rows in VCF file.
* `indSubpopIndex` index of that individual within the subpopulation
* `subpop` subpopulation ID
* `phen_sal` individual salinity phenotype
* `phen_temp` individual temperature phenotype
* `sal_opt` salinity optimum for that subpopulation
* `temp_opt` temperature optimum for that subpopulation
* `fitness` fitness of the individual in it's subpopulation
    
### 2549986039929_info.txt
Simulation information

* `seed`  simulation seed
* `sim_type` the name of the simulation file used to produce the results
* `MIG_x` migration across longitude populations
* `MIG_y` migration across latitude populations
* `N` population size
* `mu` mutation rate
* `r` recombination rate
* `BURNIN_1` length of first burnin (typically 0 opimums across metapopulation)
* `BURNIN_2` length of second burnin (typically a transition from 0 opimums across metapopulation to final environment values)
* `TOTAL_GEN` the total number of generations
* `METAPOP_SIDE_x` the number of metapopulations across longitude
* `METAPOP_SIDE_y` the number of metapopulations across latitude
* `S` Not used in this simulation.
* `SIGMA_K` Strength of stabilizing selection
* `SIGMA_QTN` Variation in QTN mutation rate

### 2549986039929_LA.txt
This file tracks the amount of local adaptation in the simulation through time

* `seed` simulation seed
* `gen` generation
* `sympatric`  average fitness of demes in sympatry
* `allopatric` average fitness of demes in allopatry
* `local_adapt` amount of local adaptation (`allopatric` - `sympatric` )
* `mean_phen0` mean of first phenotype (salinity in this simulation) across the whole metapopulation (should be near 0)
* `mean_phen1` mean of second phenotype (temperature in this simulation) across the whole metapopulation (should be near 0)
* `cor_sal_popmean` correlation between salinity and population mean phenotype (often inflated due to central limit theorem)
* `cor_temp_popmean` correlation between temperature and population mean phenotype (often inflated due to central limit theorem)
* `cor_sal_ind` correlation between salinity and individual phenotypes
* `cor_temp_ind` correlation between temperature and individual phenotypes

### 2549986039929_mig_mat.txt
This is the migration matrix that was used in the sims. It is used in pyslim for recapitation. Without it, pyslim will run forever and never recapitate.

### 2549986039929_muts.txt
This is a table of information about each mutation that is simulated SLIM AND has an allele frequency > 0.01 or < 0.99

* `seed`  simulation seed
* `mutID` mutation ID - should match mutID in other tables and VCF file
* `muttype` mutation type in SLiM
* `p` allele frequency of derived mutation - not actually the minor allele frequency
* `cor_sal` correlation of mutation allele frequency and salinity at end of simulation (all individuals) - each deme is a datapoint
* `cor_temp` correlation of mutation allele frequency and temperature at end of simulation (all individuals) - each deme is a datapoint
* `mutSalEffect` Effect of mutation on salinity
* `mutTempEffect` Effect of mutation on temperature

### 2549986039929_plusneut_MAF01.recode2.vcf.gz

* The VCF file output after pyslim and filtered to an MAF of 0.01
* Each row is a locus simulated in SLiM or added by pyslim
* Each column is an individual (all output)
* The ID of the causal mutations is in the ALT table and should match `mutID` in other tables

### 2549986039929_VCF_causal.vcf.gz

The VCF file output 

* Each row is a locus simulated in SLiM
* Each column is an individual (all output)
* The INFO field contains additional info about the mutation in SLiM that is not retained by pyslim
* This table is redundant, but used to cross-check the pyslim outputs.

### 2549986039929.trees

The trees file output by SliM that is piped into pyslim.

### 2549986039929_muts_analysis.txt

*`mut_ID` mutation ID, if not equal to 1, a causal mutation
* `seed` simulation seed
* `VCFrow` row in VCF file
*`pos_pyslim` position as output from pyslim, which is 1 less than SliM output
*`a_freq_full` allele frequency of derived allele based on all samples
*`a_freq_subset` allele frequency of derived allele based on subset of individuals sampled according to their fitness
*`mutSalEffect` for causal mutations, effect of mutation on the salinity phenotype
*`mutTempEffect`  for causal mutations, effect of mutation on the temperature phenotype
*`INFO`  for causal mutations, information output from slim
*`af_cor_temp` correlation between allele frequency and temperature based on subset of individuals sampled according to their fitness
*`af_cor_sal` correlation between allele frequency and salinity based on subset of individuals sampled according to their fitness
*`af_cor_temp_pooled` correlation between allele frequency and temperature based on subset of individuals sampled according to their fitness, for individuals pooled by temperature instead of by population
*`af_cor_sal_pooled` correlation between allele frequency and salinity based on subset of individuals sampled according to their fitness, for individuals pooled by salinity instead of by population
*`af_slope_temp` slope between allele frequency and temperature based on subset of individuals sampled according to their fitness
*`af_slope_sal` slope between allele frequency and salinity based on subset of individuals sampled according to their fitness
* `Va_temp` VA to temperature explained by locus in metapopultion (this is a little misleading because some loci are low here, but explain a lot of VA in a specific environment), based on a subset of individuals after sampling based on their fitness
* `Va_temp_prop` Proportion of total VA to temperature explained by locus in metapopultion (this is a little misleading because some loci are low here, but explain a lot of VA in a specific environment), based on a subset of individuals after sampling based on their fitness
* `Va_sal` VA to salinity explained by locus in metapopultion (this is a little misleading because some loci are low here, but explain a lot of VA in a specific environment), based on a subset of individuals after sampling based on their fitness
* `Va_sal_prop` Proportion of total VA to salinity explained by locus in metapopultion (this is a little misleading because some loci are low here, but explain a lot of VA in a specific environment), based on a subset of individuals after sampling based on their fitness

* `LEA3.2_lfmm2_mlog10P_tempenv`
* `LEA3.2_lfmm2_mlog10P_tempenv_sig`
* `LEA3.2_lfmm2_Va_temp_prop`

* `LEA3.2_lfmm2_mlog10P_salenv`
* `LEA3.2_lfmm2_mlog10P_salenv_sig`
* `LEA3.2_lfmm2_Va_sal_prop`

### 2549986039929_ind_subset_analysis.txt

* `RDA_allloci_temp_pred`
* `RDA_allloci_sal_pred`

TO DO

### 2549986039929_af_deme.txt
a matrix in which each row is a locus and each column is a deme, and the entry is the allele frequency

DO I OUTPUT THIS?

### 2549986039929_RDA_loadings.txt
RDA loadings for each environment

### 2549986039929_sim_analysis.txt

After each simulation is analyzed in R, 

*`seed` 
* `nind_samp` number of individuals in sample
* `K` number of populations used in analyses
* `all_corr_phen_temp` for all individuals, correlation between individual temp phenotype and environment temperature
* `subsamp_corr_phen_temp` after sampling 20 individuals from each patch with a probability based on their fitness, correlation individual temp phenotype and environment temperature
* `all_corr_phen_sal` for all individuals, correlation between individual sal phenotype and environment salinity
* `subsamp_corr_phen_sal` after sampling 20 individuals from each patch with a probability based on their fitness, correlation between individual sal phenotype and environment salinity
*`num_causal` number of causal loci in sim
*`num_neut` number of neutral loci in sim



* `LEA3.2_lfmm2_Va_sal_prop`  proportion of Va in salinity trait explained by LFMM outliers
* `LEA3.2_lfmm2_Va_temp_prop` proportion of Va in temperature trait explained by LFMM outliers
*`LEA3.2_lfmm2_TPR_temp` true positive rate - proportion of causal loci that are significant by LFMM
*`LEA3.2_lfmm2_TPR_sal` true positive rate - proportion of causal loci that are significant by LFMM
*`LEA3.2_lfmm2_FDR_temp` false discovery rate - proportion of outliers that are true positives
*`LEA3.2_lfmm2_FDR_sal` false discovery rate - proportion of outliers that are true positives


* `RDA_Va_temp_prop` proportion of Va in temperature trait explained by RDA outliers
* `RDA_Va_sal_prop` proportion of Va in salinity trait explained by RDA outliers
*`RDA_TPR`  true positive rate - proportion of causal loci that are significant by RDA
* `RDA_FDR` false discovery rate - proportion of outliers that are true positives



*`cor_RDA500_RDloadings_tempPhen` For an RDA based on 500 random loci, the correlation between (a linear prediction of the weighted RDA loadings) and (the individual's temp phenotype)
*`cor_RDA500_RDloadings_salPhen` For an RDA based on 500 random loci, the correlation between (a linear prediction of the weighted RDA loadings) and (the individual's salinity phenotype)

*`num_causal_sig_temp_corr` number of causal loci that have significant Spearman's correlations with temperature after Bonferroni correction
*`num_causal_sig_sal_corr` number of causal loci that have significant Spearman's correlations with salinity after Bonferroni correction
*`num_neut_sig_temp_corr` number of neutral loci that have significant  Spearman's correlations with temperature after Bonferroni correction
*`num_neut_sig_sal_corr` number of neutral loci that have significant Spearman's correlations with salinity after Bonferroni correction

*`median_causal_temp_cor` median abs(Spearman's correlation) between allele frequency and temperature for causal loci
*`median_causal_sal_cor` median abs(Spearman's correlation)  between allele frequency and salinity for causal loci
*`median_neut_temp_cor` median abs(Spearman's correlation)  between allele frequency and temperature for neutral loci
*`median_neut_sal_cor` median abs(Spearman's correlation)  between allele frequency and salinity for neutral loci
