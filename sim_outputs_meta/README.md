Below is a description of metadata for each set of output tables. The first number is the seed for the simulation.

### seed_fitnessmat_ind.txt
* an n_deme x n_individual table that indicates the fitness of each individual (in columns) in every metapopulation deme (in rows)

### seed_fitnessmat.txt
* an n_deme x n_deme table that indicates the mean fitness of individuals from the source deme (in columns) in the transplant deme (in rows) CHECK THIS IS CORRECT

### seed_ind.txt
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
    
### seed_info.txt
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

### seed_LA.txt
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

### seed_mig_mat.txt
This is the migration matrix that was used in the sims. It is used in pyslim for recapitation. Without it, pyslim will run forever and never recapitate.

### seed_muts.txt
This is a table of information about each mutation that is simulated SLIM AND has an allele frequency > 0.01 or < 0.99

* `seed`  simulation seed
* `mutID` mutation ID - should match mutID in other tables and VCF file
* `muttype` mutation type in SLiM
* `p` allele frequency of derived mutation - not actually the minor allele frequency
* `cor_sal` correlation of mutation allele frequency and salinity at end of simulation (all individuals) - each deme is a datapoint
* `cor_temp` correlation of mutation allele frequency and temperature at end of simulation (all individuals) - each deme is a datapoint
* `mutSalEffect` Effect of mutation on salinity
* `mutTempEffect` Effect of mutation on temperature

### seed_plusneut_MAF01.recode2.vcf.gz

* The VCF file output after pyslim and filtered to an MAF of 0.01 across the entire population
* Each row is a locus simulated in SLiM or added by pyslim
* Each column is an individual (all output)
* The ID of the causal mutations is in the ALT table and should match `mutID` in other tables

### seed_VCF_causal.vcf.gz

The VCF file output 

* Each row is a locus simulated in SLiM
* Each column is an individual (all output)
* The INFO field contains additional info about the mutation in SLiM that is not retained by pyslim
* This table is redundant, but used to cross-check the pyslim outputs.

### seed.trees

The trees file output by SliM that is piped into pyslim.

### seed_Rout_RDA_predictions.txt
* `nloci` the number of loci randomly sampled from the genome
* `Random_cor_temppredict_tempphen` the correlation between the RDA temperature prediction for an individual based on a random set of loci (see equation below) and the temperature phenotype of an individual
    * `temp_pred <- ind.sc[,1]*eigenvals(rdaout_rand)[1]*summary(rdaout_rand)$biplot[2,1] + ind.sc[,2]*eigenvals(rdaout_rand)[2]*summary(rdaout_rand)$biplot[2,2]`
* `Random_cor_salpredict_salphen` the correlation between the RDA salinity prediction for an individual based on a random set of loci (see equation below) and the salinity phenotype of an individual
    * `sal_pred <- ind.sc[,1]*eigenvals(rdaout_rand)[1]*summary(rdaout_rand)$biplot[1,1] + ind.sc[,2]*eigenvals(rdaout_rand)[2]*summary(rdaout_rand)$biplot[1,2]`

### seed_Rout_af_pop.txt
* Allele frequency as a function of demeID for each mutation in ADD FILE
* rows are demes
* columns are mutations in ADD FILE

### seed_Rout_af_sal.txt
* Allele frequency as a function of salinity for each mutation 
* rows are salinity values
* columns are mutations in ADD FILE

### seed_Rout_af_temp.txt
* Allele frequency as a function of temperature for each mutation
* rows are temperature values
* columns are mutations in ADD FILE

### seed_Rout_muts_full.txt TO DO

* `mut_ID` mutation ID, if not equal to 1, a causal mutation
* `seed` simulation seed
* `VCFrow` row in VCF file
* `pos_pyslim` position as output from pyslim, which is 1 less than SliM output
* `a_freq_full` allele frequency of derived allele based on all samples
* `a_freq_subset` allele frequency of derived allele based on subset of individuals sampled according to their fitness
* `mutSalEffect` for causal mutations, effect of mutation on the salinity phenotype
* `mutTempEffect`  for causal mutations, effect of mutation on the temperature phenotype
* `INFO`  for causal mutations, information output from slim
* `af_cor_temp` correlation between allele frequency and temperature based on subset of individuals sampled according to their fitness
* `af_cor_sal` correlation between allele frequency and salinity based on subset of individuals sampled according to their fitness
* `af_cor_temp_pooled` correlation between allele frequency and temperature based on subset of individuals sampled according to their fitness, for individuals pooled by temperature instead of by population
*`af_cor_sal_pooled` correlation between allele frequency and salinity based on subset of individuals sampled according to their fitness, for individuals pooled by salinity instead of by population
* `af_slope_temp` slope between allele frequency and temperature based on subset of individuals sampled according to their fitness
* `af_slope_sal` slope between allele frequency and salinity based on subset of individuals sampled according to their fitness
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

mutID                            "1"             
seed                             "1231222"       
VCFrow                           "1"             
pos_pyslim                       "5"             
a_freq_full                      "0.06025"       
a_freq_subset                    "0.0675"        
muttype                          NA              
p                                NA              
cor_sal                          NA              
cor_temp                         NA              
mutSalEffect                     NA              
mutTempEffect                    NA              
causal                           "FALSE"         
af_cor_temp                      "-0.1795185"    
af_slope_temp                    "-0.04622727"   
af_cor_temp_P                    "0.02483159"    
af_cor_sal                       "0.547292"      
af_slope_sal                     "0.1075909"     
af_cor_sal_P                     "7.850828e-12"  
af_cor_temp_mlog10P              "1.604996"      
af_cor_sal_mlog10P               "11.10508"      
causal_temp                      "neutral-linked"
causal_sal                       "neutral-linked"
LG                               "1"             
colors                           "#FFC1251A"     
Va_temp                          NA              
Va_temp_prop                     "0"             
Va_sal                           NA              
Va_sal_prop                      "0"             
cor_temp_sig                     "FALSE"         
cor_sal_sig                      "TRUE"          
LEA3.2_lfmm2_mlog10P_tempenv     "0.3529028"     
LEA3.2_lfmm2_mlog10P_tempenv_sig "FALSE"         
LEA3.2_lfmm2_mlog10P_salenv      "1.164015"      
LEA3.2_lfmm2_mlog10P_salenv_sig  "FALSE"         
structure_cor_G_LFMM_U1_modsal   "-0.1160243"    
structure_cor_G_LFMM_U1_modtemp  "-0.2870202"    
structure_cor_G_PC1              "-0.2018207"    
RDA1_score                       "0.06221698"    
RDA2_score                       "-0.07683058"   
RDA_mlog10P                      "0.7845996"     
RDA_mlog10P_sig                  "FALSE"         
af_cor_temp_pooled               "-0.6981386"    
af_cor_sal_pooled                "0.8806303"  

### seed_Rout_ind_subset.txt TO DO


seed           "1231222"   
subpopID       "1"         
indID          "5"         
indSubpopIndex "5"         
subpop         "1"         
phen_sal       "-0.928137" 
phen_temp      "-0.767099" 
sal_opt        "-1"        
temp_opt       "-1"        
fitness        "0.942323"  
subset         "TRUE"      
N              "100"       
opt0           "-1"        
opt1           "-1"        
x              "1"         
y              "1"         
PC1            "52.08"     
PC2            "12.4889"   
PC3            "39.5559"   
LFMM_U1_temp   "19.98274"  
LFMM_U1_sal    "-20.27494" 
LFMM_U2_temp   "19.94116"  
LFMM_U2_sal    "-19.56712" 
RDA1           "-2.730669" 
RDA2           "-0.8895302"
RDA_predict_tempPhen_20KSNPs
RDA_predict_salPhen_20KSNPs


### seed_Rout_simSummary.txt TO DO

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

seed                             "1231222"     
n_samp_tot                       "1000"        
n_samp_per_pop                   "10"      
sd_fitness_among_inds,
                       sd_fitness_among_pops,
                       final_LA 
K                                "3"           
Bonf_alpha                       "2.366528e-06"
numCausalLowMAFsample            "23"          
all_corr_phen_temp               "0.9604135"   
subsamp_corr_phen_temp           "0.8649487"   
all_corr_phen_sal                "0.9552496"   
subsamp_corr_phen_sal            "0.8569631"   
num_causal                       "474"         
num_non_causal                   "20654"       
num_neut                         "10410"       
num_causal_temp                  "473"         
num_causal_sal                   "473"         
LEA3.2_lfmm2_Va_temp_prop        "0"           
LEA3.2_lfmm2_Va_sal_prop         "0.06246894"  
LEA3.2_lfmm2_TPR_temp            "0"           
LEA3.2_lfmm2_TPR_sal             "0.00422833"  
LEA3.2_lfmm2_FDR_allSNPs_temp    NA            
LEA3.2_lfmm2_FDR_allSNPs_sal     "0.7142857"   
LEA3.2_lfmm2_FDR_neutSNPs_temp   NA            
LEA3.2_lfmm2_FDR_neutSNPs_sal    "0"           
LEA3.2_lfmm2_AUCPR_temp_allSNPs  "0.0281779"   
LEA3.2_lfmm2_AUCPR_temp_neutSNPs "0.06147373"  
LEA3.2_lfmm2_AUCPR_sal_allSNPs   "0.0281779"   
LEA3.2_lfmm2_AUCPR_sal_neutSNPs  "0.06147373"  
RDA_Va_temp_prop                 "0.1251184"   
RDA_Va_sal_prop                  "0.1311002"   
RDA_TPR                          "0.01476793"  
RDA_FDR_allSNPs                  "0.9688889"   
RDA_FDR_neutSNPs                 "0.9156627"   
RDA_AUCPR_allSNPs                "0.01652866"  
RDA_AUCPR_neutSNPs               "0.03505077"  
cor_RDA20000_RDloadings_tempPhen "0.8271191"   
cor_RDA20000_RDloadings_salPhen  "0.8205405"   
cor_VA_temp_prop                 "0.5516298"   
cor_VA_sal_prop                  "0.5384417"   
cor_TPR_temp                     "0.2114165"   
cor_TPR_sal                      "0.2029598"   
cor_FDR_allSNPs_temp             "0.9743655"   
cor_FDR_neutSNPs_temp            "0.9542962"   
cor_FDR_allSNPs_sal              "0.9752194"   
cor_FDR_neutSNPs_sal             "0.955535"    
cor_AUCPR_temp_allSNPs           "0.02581878"  
cor_AUCPR_temp_neutSNPs          "0.05166276"  
cor_AUCPR_sal_allSNPs            "0.02556545"  
cor_AUCPR_sal_neutSNPs           "0.05184973"  
median_causal_temp_cor           "0.2746745"   
median_causal_sal_cor            "0.2685987"   
median_neut_temp_cor             "0.217133"    
median_neut_sal_cor              "0.2153212"   
cor_PC1_temp                     "-0.5262854"  
cor_PC1_sal                      "-0.8366671"  
cor_PC2_temp                     "-0.8367879"  
cor_PC2_sal                      "0.5274318"   
cor_LFMMU1_temp                  "-0.0284181"  
cor_LFMMU1_sal                   "0.02959235"  
cor_LFMMU2_temp                  "0.02630689"  
cor_LFMMU2_sal                   "-0.04222666" 
cor_PC1_LFMMU1_temp              "0.8610124"   
cor_PC1_LFMMU1_sal               "-0.5587573"  
cor_PC2_LFMMU1_temp              "-0.5074529"  
cor_PC2_LFMMU1_sal               "-0.8283036" 



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
