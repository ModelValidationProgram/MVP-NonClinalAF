# 20210312 Initiate Repo

* transfered 20200224_Multi_m0.05_CG_2phen_2envi_5x10_coastline_keep.slim from repo MVP_prelim
and renamed it Case1_CoreSimEstuary2TraitsWithPleiotropy.slim



## Part 1:
[x] run this sim
	* added seed info to output tables
	* added simulation name to "info table"
	* added tree sequencing
	* downloaded newest SLiM version
	
* ISSUE SEED NOT PRINTING
	* if I use `paste(c(MY_SEED, MIG_x))`, then the seed prints as scientific format
	* I now use `paste(MY_SEED, MIG_x), and the seed prints in full
	
[x] run the R Markdown on the output, make sure it works

	* while I was working on this, I noticed what I thought be an error from the code I took from 
	`1-AnalyzeSims_10x5_this_has_af_env_plots_but_not_sure_about_MRFTraining_think_error_in_af_vs_env_plots`
	(I renamed the file). I'm not really sure if it was an error - I had to rewrite the code
	to be consistent with this markdown, but just wanted to flag it. The error was possibly :
	and not using the correct row when only two of the 
	three genotypes occurred at a locus. (this would have been rare)
	* in any case I fixed it, and sanity check matches SLiM output
	
In working on the R markdown, I began to think about the best way to present results -
in terms of clines (slopes) or correlations. Sometimes there is a significant correlation,
but not much of an allele frequency cline. This is a nuance I need to clarify for the paper.

I also began to think about the best way to present the results.

there are three types of false negatives:
1) both temp and sal missed
2) temp discovered; salinity missed
3) temp missed; salinity discovered

With correlations in af these three are easy to distinguish. But with RDA it's could be bit harder 
to characterize. If each variable loads onto an RDA axis, then I can look for significance
with that axis - I'm not sure this would tell me anything new. In addition, a locus could
be an outlier by the Mahalanobis distance metric without loading onto any one axis significantly.

thought: it's not about discovery - it's about whether discovery can tell us about the ways
	alleles confer an advantage to the environment

(i) just look at distribution of allele frequency correlations?
(ii) When do correlations correctly infer the way alleles confer a fitness advantage in an environment?
	- correlation outliers (temp, sal)
	- evaluate the 3 false positives
(iii)  When does RDA correctly infer the way alleles confer a fitness advantage in an environment?
	- RDA outliers (temp, sal)
	- I wonder if there is a way to predict the environment from the redundancy axis

(ii)  if we take the loci with significant correlations with temperature, how much of 
the additive genetic variance in the trait is explained? For example, if we take the loci with significant correlations with temperature, how much of 
the additive genetic variance in the trait is explained?
	- this is straightforward for corrs, not so much for RDA

I think I need a population-specific approach - for 6 categories of populations, in each one calculate which alleles are needed to explain the VA in that pop.
- high salinity high temp
- low salinity high temp
- high salinity med temp
- low salinity med temp
- high salinity intermediate temp
- low salinity intermediate temp

## RDA
[x] run the RDA on data I have, see what it gives

Results are interesting. When I model the genetic basis perfectly, the loci scores are crap, but the individual loadings can accurately predict an individual's temp or salinity phenotype. 

[x] add sub pop xy locations to output
	

## Plan:	
[] add neutral mutations with pyslim
[] show that RDA can predict individuals environments when we know the adaptive architecture perfectly
[] ID outliers with LFMM - TEST IF RDA CAN PREDICT IND ENVIIRONMENTS
[] ID outliers with RDA - TEST IF RDA CAN PREDICT IND ENVIIRONMENTS
[] ID outliers with OUTFLANK - TEST IF RDA CAN PREDICT IND ENVIIRONMENTS

Part 2:
- Expand this case to the 2 traits without pleiotropy
	* be sure to add simulation name to "info table"
- Expand this case to the 1 trait without pleiotropy