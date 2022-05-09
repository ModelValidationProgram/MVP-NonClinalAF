
### seed_Rout_ind_subset.txt  

(output after analysis of simulations) 
A table with information about each individual that was sampled for analysis (10 ind/deme x 100 demes)

* `seed`    simulation seed
* `subpopID` ID of deme in SliM sim
* `indID`   individual ID for the subset of individuals
* `indSubpopIndex` index of individual within that subpop    
* `subpop`        redundant with `subpopID`   
* `phen_sal`     salinity (Env2) phenotype
* `phen_temp`    temperature phenotype
* `sal_opt`      optimum salinity (Env2) of the deme where it was sampled    
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
* `LFMM_U1_sal`    loading of individual on 1st latent factor from LFMM salinity (Env2) model
* `LFMM_U2_temp`    loading of individual on 2nd latent factor from LFMM temp model
* `LFMM_U2_sal`    loading of individual on 2nd latent factor from LFMM salinity (Env2) model
* `RDA1`           loading of individual on first RDA axis for RDA model: genotype~Env
* `RDA2`           loading of individual on second RDA axis for RDA model: genotype~Env
* `RDA1_corr` 	loading of individual on first RDA axis for RDA model with structure correction: genotype~Env + Condition(PC1 + PC2)
* `RDA2_corr`  loading of individual on second RDA axis for RDA model with structure correction: genotype~Env + Condition(PC1 + PC2)
* `RDA_predict_tempPhen_20KSNPs`temperature phenotype prediction from an RDA with 20K random loci. RDA model: genotype~Env. See `seed_Rout_RDA_predictions.txt`
* `RDA_predict_salPhen_20KSNPs`salinity (Env2) phenotype prediction from an RDA with 20K random loci. RDA model: genotype~Env.  See `seed_Rout_RDA_predictions.txt`
* `RDA_predict_tempPhen_20KSNPs_structcorr`temperature phenotype prediction from an RDA with structure correction and 20K random loci. RDA model: genotype~Env + Condition(PC1 + PC2). See `seed_Rout_RDA_predictions.txt`
* `RDA_predict_salPhen_20KSNPs_structcorr` salinity (Env2) phenotype prediction from an RDA with structure correction and 20K random loci. RDA model: genotype~Env+ Condition(PC1 + PC2).  See `seed_Rout_RDA_predictions.txt`
