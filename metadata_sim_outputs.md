Below is a description of metadata for each set of output tables. The first number is the seed for the simulation.

Note: The origin of the simulations was to better understand how oysters adapt to salinity and temperature gradients in estuaries. Opt0 in the simulations corresponds to the salinity trait, which also corresponds to Env2 in the publications.

### seed_fitnessmat_ind.txt (SLiM output)
* an n_deme x n_individual table that indicates the fitness of each individual (in columns) in every metapopulation deme (in rows)

### seed_fitnessmat.txt (SLiM output)
* an n_deme x n_deme table that indicates the mean fitness of individuals from the source deme (in columns) in the transplant deme (in rows) CHECK THIS IS CORRECT

### seed_ind.txt (SLiM output) Individual information output last generation

* `seed` simulation seed
* `indID` individual ID in SliM. Corresponds to rows in VCF file.
* `indSubpopIndex` index of that individual within the subpopulation
* `subpop` subpopulation ID
* `phen_sal` individual salinity phenotype
* `phen_temp` individual temperature phenotype
* `sal_opt` salinity optimum for that subpopulation
* `temp_opt` temperature optimum for that subpopulation
* `fitness` fitness of the individual in it's subpopulation
    
### seed_info.txt (SLiM output) Simulation information

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

### seed_LA.txt (SLiM output) This file tracks the amount of local adaptation in the simulation through time

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

### seed_mig_mat.txt (SLiM output)
This is the migration matrix that was used in the sims. It is used in pyslim for recapitation. Without it, pyslim will run forever and never recapitate.

### seed_muts.txt (SLiM output) mutation table
This is a table of information about each mutation that is simulated SLIM AND has an allele frequency > 0.01 or < 0.99

* `seed`  simulation seed
* `mutID` mutation ID - should match mutID in other tables and VCF file
* `muttype` mutation type in SLiM
* `p` allele frequency of derived mutation - not actually the minor allele frequency
* `cor_sal` correlation of mutation allele frequency and salinity at end of simulation (all individuals) - each deme is a datapoint
* `cor_temp` correlation of mutation allele frequency and temperature at end of simulation (all individuals) - each deme is a datapoint
* `mutSalEffect` Effect of mutation on salinity
* `mutTempEffect` Effect of mutation on temperature

### seed_popInfo.txt (SLiM output) deme information
* `seed`  simulation seed
* `subpopID`
* `N` number of individuals simmulated in that deme
* `opt0` salinity optimum
* `opt1` Temperature optimum
* `x` x location on grid
* `y` y location on grid

### seed_plusneut_MAF01.recode2.vcf.gz

* The VCF file output after pyslim and filtered to an MAF of 0.01 across the entire population
* Each row is a locus simulated in SLiM or added by pyslim
* Each column is an individual (all output)
* The ID of the causal mutations is in the ALT table and should match `mutID` in other tables

### seed_VCF_causal.vcf.gz (SLiM output) The VCF file output from SLiM

* Each row is a locus simulated in SLiM
* Each column is an individual (all output)
* The INFO field contains additional info about the mutation in SLiM that is not retained by pyslim
* This table is redundant, but used to cross-check the pyslim outputs.

### seed.trees (pyslim output)

The trees file output by SliM that is piped into pyslim.

### seed_Rout_RDA_predictions.txt (output after analysis of simulations)

correlation between (a linear prediction of the weighted RDA loadings) and (the individual's temp phenotype)

* `nloci` the number of loci randomly sampled from the genome
* `Random_cor_temppredict_tempphen` the correlation between the RDA temperature prediction for an individual based on a random set of loci (see equation below) and the temperature phenotype of an individual
    * `temp_pred <- ind.sc[,1]*eigenvals(rdaout_rand)[1]*summary(rdaout_rand)$biplot[2,1] + ind.sc[,2]*eigenvals(rdaout_rand)[2]*summary(rdaout_rand)$biplot[2,2]`
* `Random_cor_salpredict_salphen` the correlation between the RDA salinity prediction for an individual based on a random set of loci (see equation below) and the salinity phenotype of an individual
    * `sal_pred <- ind.sc[,1]*eigenvals(rdaout_rand)[1]*summary(rdaout_rand)$biplot[1,1] + ind.sc[,2]*eigenvals(rdaout_rand)[2]*summary(rdaout_rand)$biplot[1,2]`

### seed_Rout_af_pop.txt (output after analysis of simulations)
* Allele frequency as a function of demeID for each mutation in ADD FILE
* rows are demes that correspond to the popInfo.txt
* columns are mutations in ADD FILE

### seed_Rout_af_sal.txt (output after analysis of simulations)
* Allele frequency as a function of salinity for each mutation 
* rows are salinity values
* columns are mutations in ADD FILE

### seed_Rout_af_temp.txt (output after analysis of simulations)
* Allele frequency as a function of temperature for each mutation
* rows are temperature values
* columns are mutations in ADD FILE

### `seed_Rout_muts_full.txt` (output after analysis of simulations) A table with information about each locus

* `mut_ID` mutation ID, if not equal to 1, a causal mutation
* `seed` simulation seed
* `VCFrow` row in VCF file
* `pos_pyslim` position as output from pyslim, which is 1 less than SliM output
* `a_freq_full` allele frequency of derived allele based on all samples
* `a_freq_subset` allele frequency of derived allele based on subset of individuals sampled according to their fitness
* `muttype` muttype from SliM (causal mutations)
* `p` allele frenquency from SliM (causal mutations)
* `cor_sal` correlation of deme allele frequency and deme salinity optimum from SliM (all samples)
* `cor_temp` correlation of deme allele frequency and deme temperature optimum from SliM (all samples)
* `mutSalEffect` for causal mutations, effect of mutation on the salinity phenotype
* `mutTempEffect`  for causal mutations, effect of mutation on the temperature phenotype
* `causal`a logical value indiciating whether it was a causal mutation (had a non-zero effect on either temp or salinity trait)

* `INFO`  for causal mutations, information output from slim

* `af_cor_temp` 

* `af_cor_temp` correlation between allele frequency and temperature based on subset of individuals sampled according to their fitness
* `af_slope_temp` slope between between allele frequency and temperature based on subset of individuals sampled according to their fitness
* `af_cor_temp_P` P-value for correlation between allele frequency and temperature based on subset of individuals sampled according to their fitness
* `af_cor_temp_mlog10P` -log10 P value of above measure
* `af_cor_sal` correlation between allele frequency and salinity based on subset of individuals sampled according to their fitness
* `af_slope_sal` slope between between allele frequency and salinity based on subset of individuals sampled according to their fitness
* `af_cor_sal_P`  P-value for correlation between allele frequency and salinity based on subset of individuals sampled according to their fitness
* `af_cor_sal_mlog10P` -log10 P value of above measure    

* `causal_temp` a categorical variable indicating the locus type (see sub-bullets below)
* `causal_sal` a categorical variable if the locus type  (see sub-bullets below)
    * `causal`  the locus has an allele with a non-zero effect on the phenotype
    * `neutral-linked` a locus with alleles that have zero effect on the phenotype, but arise in the part of the genome where they may be linked to a locus that is causal. (e.g., sites 1-500000)
    * `neutral` a locus with alleles that have zero effect on the phenotype, and arise in the part of the genome unaffected by selection (e.g., sites 500001-1e06)
* `LG` linkage group       
* `colors` a color level for plotting - shades of yellow for LGs affected by selection(e.g., sites 1-500000), shades of grey for LGs unaffected by selection (e.g., sites 500001-1e06)

* `va_temp_full` true value additive genetic variance of the locus on the temeprature trait, based on entire sample
* `va_sal_full` true value additive genetic variance of the locus on the temeprature trait, based on entire sample
* `va_temp_full_prop` the true proportion of Va explained by this locus on the temperature trait, based on the entire sample
* `va_sal_full_prop` the true proportion of Va explained by this locus on the salinity trait, based on the entire sample
* `Va_temp` estimate of additive genetic variance of the locus on the temeprature trait (based on the subset of individuals sampled according to their fitness)     
* `Va_temp_prop` proportion of the total additive genetic variance on the temeprature trait explained by this mutation      
* `Va_sal` estimate of additive genetic variance of the locus on the salinity trait (based on the subset of individuals sampled according to their fitness)     
* `Va_sal_prop` proportion of the total additive genetic variance on the salinity trait explained by this mutation      

* `He_outflank` Expected heterozygosity based on the sample, calculated in outflank
* `Fst_outflank` W&C's 1984 Fst based on the sample, calculated in outflank
* 

* `cor_temp_sig`  a logical value indicating whether the `af_cor_temp` was significant after Bonferroni correction       
* `cor_sal_sig`  a logical value indicating whether the `af_cor_sal` was significant after Bonferroni correction 

* `af_cor_temp_pooled` correlation between allele frequency and temperature based on subset of individuals sampled according to their fitness, for individuals pooled by temperature instead of by population
*`af_cor_sal_pooled` correlation between allele frequency and salinity based on subset of individuals sampled according to their fitness, for individuals pooled by salinity instead of by population
* `af_slope_temp` slope between allele frequency and temperature based on subset of individuals sampled according to their fitness
* `af_slope_sal` slope between allele frequency and salinity based on subset of individuals sampled according to their fitness
* `Va_temp` VA to temperature explained by locus in metapopultion (this is a little misleading because some loci are low here, but explain a lot of VA in a specific environment), based on a subset of individuals after sampling based on their fitness
* `Va_temp_prop` Proportion of total VA to temperature explained by locus in metapopultion (this is a little misleading because some loci are low here, but explain a lot of VA in a specific environment), based on a subset of individuals after sampling based on their fitness
* `Va_sal` VA to salinity explained by locus in metapopultion (this is a little misleading because some loci are low here, but explain a lot of VA in a specific environment), based on a subset of individuals after sampling based on their fitness
* `Va_sal_prop` Proportion of total VA to salinity explained by locus in metapopultion (this is a little misleading because some loci are low here, but explain a lot of VA in a specific environment), based on a subset of individuals after sampling based on their fitness

* `LEA3.2_lfmm2_mlog10P_tempenv` -log10 P value from the lfmm model for temperature
* `LEA3.2_lfmm2_mlog10P_tempenv_sig` a logical value indicating if the q-value from the lfmm model for temperature was less than 0.05
* `LEA3.2_lfmm2_mlog10P_salenv` -log10 P value from the lfmm model for salinity
* `LEA3.2_lfmm2_mlog10P_salenv_sig` a logical value indicating if the q-value from the lfmm model for salinity was less than 0.05
* `structure_cor_G_LFMM_U1_modsal` for this locus, correlation between genotypes and the individual loadings on the first latent factor of the lfmm salinity model. (e.g., is this locus correlated with what the model thinks is structure)
* `structure_cor_G_LFMM_U1_modtemp` for this locus, correlation between genotypes and the individual loadings on the first latent factor of the lfmm temperature model. (e.g., is this locus correlated with what the model thinks is structure)

* `structure_cor_G_PC1` for this locus, correlation between genotypes and the individual loadings on the first PC axis. (e.g., is this locus correlated with actual structure)

* `RDA1_score` loading of the locus on the first RDA axis
* `RDA2_score` loading of the locus on the second RDA axis
* `RDA_mut_sal_pred` RDA prediction of effect mutation has on salinity
* `RDA_mut_temp_pred` RDA prediction of effect mutation has on temperature
* `RDA_mlog10P` -log 10 P-value from the "rdadapt" function from Capblanq et al. 20XX
* `RDA_mlog10P_sig` a logical variable indicating if the RDA q-value was less than 0.001 (a more conservative estimate was used than with lfmm because lack of structure correction led to a large number of false positive results)
* `af_cor_temp_pooled` correlation between allele frequency and temperature if all demes were pooled together according to their temperature (this inflates correlations)  
* `af_cor_sal_pooled` correlation between allele frequency and salinity if all demes were pooled together according to their salinity (this inflates correlations)

### seed_Rout_ind_subset.txt  (output after analysis of simulations) A table with information about each individual that was sampled for analysis


* `seed`    simulation seed
* `subpopID` ID of deme in SliM sim
* `indID`   individual ID for the subset of individuals
* `indSubpopIndex` index of individual within that subpop    
* `subpop`        redundant with `subpopID`   
* `phen_sal`     salinity phenotype
* `phen_temp`    temperature phenotype
* `sal_opt`      optimum salinity of the deme where it was sampled    
* `temp_opt`      optimum temperature of the deme where it was sampled    
* `fitness`       fitness in the deme where it was sampled
* `subset`       logical indicating if the individual was included in all analyses (should be all TRUE)    
* `N`             population size of the deme that the individual was sampled from      
* `opt0`          redundant with `sal_opt`       
* `opt1`           redundant with `temp_opt`     
* `x`             x location of the deme where it was sampled     
* `y`             y location of the deme where it was sampled    
* `PC1`          loading of individual on PC1 axis    
* `PC2`             loading of individual on PC2 axis 
* `PC3`        loading of individual on PC3 axis
* `LFMM_U1_temp`    loading of individual on 1st latent factor from LFMM temp model 
* `LFMM_U1_sal`    loading of individual on 1st latent factor from LFMM salinity model
* `LFMM_U2_temp`    loading of individual on 2nd latent factor from LFMM temp model
* `LFMM_U2_sal`    loading of individual on 2nd latent factor from LFMM salinity model
* `RDA1`           loading of individual on first RDA axis
* `RDA2`           loading of individual on second RDA axis
* `RDA_predict_tempPhen_20KSNPs`temperature phenotype prediction from an RDA with 20K random loci. See `seed_Rout_RDA_predictions.txt`
* `RDA_predict_salPhen_20KSNPs`salinity phenotype prediction from an RDA with 20K random loci.  See `seed_Rout_RDA_predictions.txt`


### `seed_Rout_simSummary.txt` (output after analysis of simulations) A table with information about the entire simulation

After each simulation is analyzed in R, 

*`seed` simulation seed
* `n_samp_tot` total number of individuals sampled
* `nind_samp_per_pop` number of individuals sampled from each deme
* `sd_fitness_among_inds` variance in fitness among all sampled individuals in the simulation (sampling prob. is proportional to fitness to mimic viability selection)
* `sd_fitness_among_pops` variance in fitness among all demes in the simulation after sampling (sampling prob. is proportional to fitness to mimic viability selection)
* `final_LA` final amount of local adaptation in the simulation
* `K` number of populations used in analyses
* `Bonf_alpha` the significance threshold for P-values applied to the correlation
* `numCausalLowMAFsample` number of causal loci that were not filtered out, but were below the MAF cutoff. These were included in the calculations. **REVISIT**
* `all_corr_phen_temp` for all individuals, correlation between individual temp phenotype and environment temperature
* `subsamp_corr_phen_temp` after sampling 20 individuals from each patch with a probability based on their fitness, correlation individual temp phenotype and environment temperature
* `all_corr_phen_sal` for all individuals, correlation between individual sal phenotype and environment salinity
* `subsamp_corr_phen_sal` after sampling 20 individuals from each patch with a probability based on their fitness, correlation between individual sal phenotype and environment salinity
* `num_causal` number of causal loci in sim
* `num_non_causal`number of neutral loci in sim arising on the half of the genome where they could be linked to causal loci
* `num_neut` number of truly neutral loci in sim
* `num_causal_temp` number of loci with non-zero phenotypic effects on the temperature phenotype
* `num_causal_sal`number of loci with non-zero phenotypic effects on the salinity phenotype

* `va_temp_total` total additive genetic variance in the temperatuer trait, based on the entire sample
* `va_sal_total` total additive genetic variance in the salinity trait, based on the entire sample
* `mean_Fst` overall FST calculated from mean(T1)/mean(T2) in outflank
----

* `LEA3.2_lfmm2_Va_temp_prop` proportion of additive genetic variance (Va) in the temperature trait explained by outliers in the LFMM temp model
* `LEA3.2_lfmm2_Va_sal_prop` proportion of additive genetic variance (Va) in the saliniity trait explained by outliers in the LFMM salinity model
*  `LEA3.2_lfmm2_TPR_temp` true positive rate of the LFMM temp model for loci with alleles that have non-zero effects on the temperature phenotype
*  `LEA3.2_lfmm2_TPR_sal` true positive rate of the LFMM salinity model for loci with alleles that have non-zero effects on the salnity phenotype
*  `LEA3.2_lfmm2_FDR_allSNPs_temp` false discovery rate of the LFMM temp model for the entire genome      
* `LEA3.2_lfmm2_FDR_allSNPs_sal` false discovery rate of the LFMM temp model for the entire genome   
* `LEA3.2_lfmm2_FDR_neutSNPs_temp` an optimistic calculation of the false discovery rate of the LFMM temp model, including only causal loci and neutral loci not affected by selection (any non-causal loci that arises on the half of the genome affected by selection was excluded)
* `LEA3.2_lfmm2_FDR_neutSNPs_sal` an optimistic calculation of the false discovery rate of the LFMM salinity model, including only causal loci and neutral loci not affected by selection (any non-causal loci that arises on the half of the genome affected by selection was excluded)
*  `LEA3.2_lfmm2_AUCPR_temp_allSNPs` the AUC-PR of the lfmm temp model based on the entire genome  
* `LEA3.2_lfmm2_AUCPR_temp_neutSNPs` an optimistic calculation of the AUC-PR of the LFMM temp model, including only causal loci and neutral loci not affected by selection (any non-causal loci that arises on the half of the genome affected by selection was excluded)
* `LEA3.2_lfmm2_AUCPR_sal_allSNPs` the AUC-PR of the lfmm temp model based on the entire genome
* `LEA3.2_lfmm2_AUCPR_sal_neutSNPs` the AUC-PR of the lfmm salinity model based on the entire genome  

----

* `RDA_Va_temp_prop` proportion of additive genetic variance (Va) in the temperature trait explained by outliers in the RDA outlier analysis, following (Capblanq 2018, https://onlinelibrary.wiley.com/doi/abs/10.1111/1755-0998.12906)
* `RDA_Va_sal_prop` proportion of additive genetic variance (Va) in the salinity trait explained by outliers in the RDA outlier analysis, following (Capblanq 2018, https://onlinelibrary.wiley.com/doi/abs/10.1111/1755-0998.12906)
* `RDA_TPR` true positive rate of the RDA for all causal loci. Since RDA is a multidimensional analysis, I did not differentiate between loci that had causal effects on temperature or salinity.
* `RDA_FDR_allSNPs` false discovery rate of the RDA outlier analysis based on the entire genome
* `RDA_FDR_neutSNPs` an optimistic calculation of the false discovery rate of the RDA, including only causal loci and neutral loci not affected by selection (any non-causal loci that arises on the half of the genome affected by selection was excluded)
* `RDA_AUCPR_allSNPs`AUC-PR of the RDA outlier analysis based on the entire genome
* `RDA_AUCPR_neutSNPs` an optimistic calculation of the AUC-PR of the RDA, including only causal loci and neutral loci not affected by selection (any non-causal loci that arises on the half of the genome affected by selection was excluded)
* `cor_RDA20000_RDloadings_tempPhen` correlation between the true temperature phenotype and that predicted from an RDA based on 20K SNPs. see `seed_Rout_RDA_predictions`
* `cor_RDA20000_RDloadings_salPhen` correlation between the true salinity phenotype and that predicted from an RDA based on 20K SNPs. see `seed_Rout_RDA_predictions`  

----

* `cor_VA_temp_prop` proportion of VA in temperature phenotype explained by clinal outliers for temperature, based on kendall's correlation between deme allele frequency and deme temperature
* `cor_VA_sal_prop`  proportion of VA in salinity phenotype explained by clinal outliers for salinity, based on kendall's correlation between deme allele frequency and deme salinity
* `cor_TPR_temp` true positive rate for loci with non-zero effects on temperature, based on kendall's correlation between deme allele frequency and deme temperature 
* `cor_TPR_sal`true positive rate for loci with non-zero effects on salinity, based on kendall's correlation between deme allele frequency and deme salinity
* `cor_FDR_allSNPs_temp` false discovery rate of (kendall's correlation between deme allele frequency and deme temperature) for loci with non-zero effects on temperature 
* `cor_FDR_neutSNPs_temp` an optimistic calculation for false discovery rate of  (kendall's correlation between deme allele frequency and deme temperature) for loci with non-zero effects on temperature, excluding non-causal loci in half of genome affected by selection 
* `cor_FDR_allSNPs_sal` false discovery rate of (kendall's correlation between deme allele frequency and deme salinity) for loci with non-zero effects on salinity
* `cor_FDR_neutSNPs_sal` an optimistic calculation for false discovery rate of  (kendall's correlation between deme allele frequency and deme temperature) for loci with non-zero effects on temperature, excluding non-causal loci in half of genome affected by selection
* `cor_AUCPR_temp_allSNPs` AUC-PR of kendall's correlation between deme allele frequency and deme temperature, based on the whole genome and causal loci for temperature
* `cor_AUCPR_temp_neutSNPs` an optimistic estimate of AUC-PR of kendall's correlation between deme allele frequency and deme temperature, based on causal loci for temperature and neutral loci not affected by selection (excluding non-causal loci in half of genome affected by selection)
* `cor_AUCPR_sal_allSNPs` AUC-PR of kendall's correlation between deme allele frequency and deme salinity, based on the whole genome and causal loci for salinity
* `cor_AUCPR_sal_neutSNPs` an optimistic estimate of AUC-PR of kendall's correlation between deme allele frequency and deme salinity, based on causal loci for salinity and neutral loci not affected by selection (excluding non-causal loci in half of genome affected by selection)

----

* `median_causal_temp_cor` median abs(Spearman's correlation) between allele frequency and temperature for causal loci
* `median_causal_sal_cor` median abs(Spearman's correlation)  between allele frequency and salinity for causal loci
* `median_neut_temp_cor` median abs(Spearman's correlation)  between allele frequency and temperature for neutral loci
* `median_neut_sal_cor` median abs(Spearman's correlation)  between allele frequency and salinity for neutral loci

----

* `cor_PC1_temp` correlation between individual loading on PC1 from the principle components based on the Genotype-matrix (individual genotypes labeled as 0,1,2) and temperature of the deme where it was sampled
* `cor_PC1_sal` correlation between individual loading on PC1 from the principle components based on the Genotype-matrix and salnity of the deme where it was sampled
* `cor_PC2_temp` correlation between individual  loading on  PC2 from the principle components based on the Genotype-matrix and temperature of the deme where it was sampled
* `cor_PC2_sal` correlation between individual  loading on PC2 from the principle components based on the Genotype-matrix and salnity of the deme where it was sampled
* `cor_LFMMU1_temp` correlation between the individual loading on the latent factor 1 from the lfmm model based on temperature
* `cor_LFMMU1_sal` correlation between the individual loading on the latent factor 1 from the lfmm model based on salnity
* `cor_LFMMU2_temp` correlation between the individual loading on the latent factor 2 from the lfmm model based on temperature
* `cor_LFMMU2_sal` correlation between the individual loading on the latent factor 2 from the lfmm model based on salnity
* `cor_PC1_LFMMU1_temp` correlation between (individual loading on PC1 from the principle components based on the Genotype-matrix) and (individual loading on the latent factor 1 from the lfmm model based on temperature)
* `cor_PC1_LFMMU1_sal` correlation between (individual loading on PC1 from the principle components based on the Genotype-matrix) and (individual loading on the latent factor 1 from the lfmm model based on salinity)
* `cor_PC2_LFMMU1_temp` correlation between (individual loading on PC2 from the principle components based on the Genotype-matrix) and (individual loading on the latent factor 1 from the lfmm model based on temperature)
* `cor_PC2_LFMMU1_sal` correlation between (individual loading on PC2 from the principle components based on the Genotype-matrix) and (individual loading on the latent factor 1 from the lfmm model based on salinity)

----

* `cor_af_temp_noutliers` number of outliers for cor(af,temp) after Bonferroni correction
* `cor_af_sal_noutliers` number of outliers for cor(af,salinity) after Bonferroni correction
* `nSNPs` total number of SNPs in analysis
* `cor_FPR_temp_neutSNPs` false positive rate in cor(af,temp) after Bonferroni correction, based on loci unaffected by selection
* `cor_FPR_sal_neutSNPs` false positive rate in cor(af,sal) after Bonferroni correction, based on loci unaffected by selection
* `LEA3.2_lfmm2_mlog10P_tempenv_noutliers` number of outliers for the lfmm temp model (qvalue <0.05)
* `LEA3.2_lfmm2_mlog10P_salenv_noutliers` number of outliers for the lfmm salinity model (qvalue <0.05)
* `LEA3.2_lfmm2_num_causal_sig_temp` number of causal loci on the temp trait, significant in the lfmm temp model (qvalue <0.05)
* `LEA3.2_lfmm2_num_neut_sig_temp` number of neutral loci false positives (only neutral loci not affected by selection), significant in the lfmm temp model (qvalue <0.05)
* `LEA3.2_lfmm2_num_causal_sig_sal` number of causal loci on the salinity trait, significant in the lfmm salinity model (qvalue <0.05)
* `LEA3.2_lfmm2_num_neut_sig_sal` number of neutral loci false positives (only neutral loci not affected by selection), significant in the lfmm salinity model (qvalue <0.05)
*  `LEA3.2_lfmm2_FPR_neutSNPs_temp` false positive rate of lfmm temperature model
*  `LEA3.2_lfmm2_FPR_neutSNPs_sal` false positive rate of lfmm salinity model

----

*  `RDA1_propvar` proportion of variance explained by first RDA axis
*  `RDA2_propvar` proportion of variance explained by second RDA axis
*  `RDA1_temp_cor` output of `summary(rdaout)$biplot[2,1]`, which is the correlation between RDA1 and the temperature environmental variable
*  `RDA1_sal_cor`  output of `summary(rdaout)$biplot[1,1]`, which is the correlation between RDA1 and the salinity environmental variable
*  `RDA2_temp_cor` output of `summary(rdaout)$biplot[2,2]`, which is the correlation between RDA2 and the temperature environmental variable
*  `RDA2_sal_cor` output of `summary(rdaout)$biplot[1,2]`, which is the correlation between RDA2 and the salinity environmental variable
*  `RDA_mlog10P_sig_noutliers` number of outliers in the RDA analysis 
*  `RDA_FPR_neutSNPs` false positive rate of the RDA analysis, based only on neutral SNPs
*  `RDA_RDALoading_cor_tempEffect` correlation between RDA prediction of a mutation effect on temperature and the true effect size of allele on temperature
*  `RDA_RDALoading_cor_salEffect` correlation between RDA prediction a mutation effect on salinity and the true effect size of allele on salinity


