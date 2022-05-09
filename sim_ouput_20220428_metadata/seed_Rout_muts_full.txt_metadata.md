
### `seed_Rout_muts_full.txt` 

A table with information about each locus used for analysis (loci after sampling individuals
and then filtering MAF > 0.01)

* `mut_ID` mutation ID, if not equal to 1, a causal mutation
* `seed` simulation seed
* `VCFrow` row in VCF file
* `LG` linkage group
* `pos_pyslim` position as output from pyslim, which is 1 less than SliM output
* `mutname` A unique ID for the mutation based on linkage group and position. Note that in SLiM,
a few rare mutations may arise on the same "position" on the SLiM genetic map but at differen times.
In analysis, these are treated as different loci that are so closely linked that they are below the
resolution of the SLiM genetic map.
* `a_freq_full` allele frequency of derived allele based on all samples
* `a_freq_subset` allele frequency of derived allele based on subset of individuals sampled according to their fitness
* `muttype` muttype from SliM (causal mutations)
* `p` allele frenquency from SliM (causal mutations)
* `cor_sal` correlation of deme allele frequency and deme salinity optimum from SliM (all samples)
* `cor_temp` correlation of deme allele frequency and deme temperature optimum from SliM (all samples)
* `mutSalEffect` for causal mutations, effect of mutation on the salinity phenotype
* `mutTempEffect`  for causal mutations, effect of mutation on the temperature phenotype

* `va_temp_full` true value additive genetic variance of the locus on the temeprature trait, based on entire sample
* `va_sal_full` true value additive genetic variance of the locus on the temeprature trait, based on entire sample
* `va_temp_full_prop` the true proportion of Va explained by this locus on the temperature trait, based on the entire sample
* `va_sal_full_prop` the true proportion of Va explained by this locus on the salinity trait, based on the entire sample

* `causal`a logical value indiciating whether it was a causal mutation (had a non-zero effect on either temp or salinity trait)

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
* `colors` a color level for plotting - shades of yellow for LGs affected by selection(e.g., sites 1-500000), shades of grey for LGs unaffected by selection (e.g., sites 500001-1e06)

* `Va_temp` estimate of additive genetic variance of the locus on the temeprature trait (based on the subset of individuals sampled according to their fitness)     
* `Va_temp_prop` proportion of the total additive genetic variance on the temeprature trait explained by this mutation      
* `Va_sal` estimate of additive genetic variance of the locus on the salinity trait (based on the subset of individuals sampled according to their fitness)     
* `Va_sal_prop` proportion of the total additive genetic variance on the salinity trait explained by this mutation      

* `cor_temp_sig`  a logical value indicating whether the `af_cor_temp` was significant after Bonferroni correction       
* `cor_sal_sig`  a logical value indicating whether the `af_cor_sal` was significant after Bonferroni correction 


* `He_outflank` Expected heterozygosity based on the sample, calculated in outflank
* `Fst_outflank` W&C's 1984 Fst based on the sample, calculated in outflank


* `LEA3.2_lfmm2_mlog10P_tempenv` -log10 P value from the lfmm model for temperature
* `LEA3.2_lfmm2_mlog10P_tempenv_sig` a logical value indicating if the q-value from the lfmm model for temperature was less than 0.05
* `LEA3.2_lfmm2_mlog10P_salenv` -log10 P value from the lfmm model for salinity
* `LEA3.2_lfmm2_mlog10P_salenv_sig` a logical value indicating if the q-value from the lfmm model for salinity was less than 0.05


* `structure_cor_G_LFMM_U1_modsal` for this locus, correlation between genotypes and the individual loadings on the first latent factor of the lfmm salinity model. (e.g., is this locus correlated with what the model thinks is structure)
* `structure_cor_G_LFMM_U1_modtemp` for this locus, correlation between genotypes and the individual loadings on the first latent factor of the lfmm temperature model. (e.g., is this locus correlated with what the model thinks is structure)

* `structure_cor_G_PC1` for this locus, correlation between genotypes and the individual loadings on the first PC axis. (e.g., is this locus correlated with actual structure)

* `RDA1_score` loading of the locus on the first RDA axis. RDA model: genotype ~ environment
* `RDA2_score` loading of the locus on the second RDA axis. RDA model: genotype ~ environment
* `RDA1_score_corr` loading of the locus on the first RDA axis. RDA model with structure correction: genotype ~ environment + Condition(PC1 + PC2)
* `RDA2_score_corr` loading of the locus on the second RDA axis. RDA model with structure correction: genotype ~ environment + Condition(PC1 + PC2)

* `RDA_mlog10P` -log 10 P-value from the `rdadapt` function from Capblanq et al. RDA model: genotype ~ environment
* `RDA_mlog10P_sig` a logical variable indicating if the RDA q-value was less than 0.05. RDA model: genotype ~ environment
* `RDA_mlog10P_corr` -log 10 P-value from the `rdadapt` function from Capblanq et al. RDA model with structure correction: genotype ~ environment + Condition(PC1 + PC2)
* `RDA_mlog10P_sig_corr` a logical variable indicating if the RDA q-value was less than 0.05. RDA model with structure correction: genotype ~ environment + Condition(PC1 + PC2)


* `RDA_mut_sal_pred` RDA prediction of effect mutation has on salinity. RDA model: genotype ~ environment
* `RDA_mut_temp_pred` RDA prediction of effect mutation has on temperature. RDA model: genotype ~ environment
* `RDA_mut_sal_pred_structcorr` RDA prediction of effect mutation has on salinity. RDA model with structure correction: genotype ~ environment + Condition(PC1 + PC2)
* `RDA_mut_temp_pred_structcorr` RDA prediction of effect mutation has on temperature. RDA model with structure correction: genotype ~ environment + Condition(PC1 + PC2)

* `af_cor_temp_pooled` correlation between allele frequency and temperature if all demes were pooled together according to their temperature (this inflates correlations)  
* `color_af.temp.cline` color based on the previous variable for plotting
* `af_cor_sal_pooled` correlation between allele frequency and salinity if all demes were pooled together according to their salinity (this inflates correlations)
* `color_af.sal.cline` color based on the previous variable for plotting


* `INFO`  for causal mutations, information output from slim
