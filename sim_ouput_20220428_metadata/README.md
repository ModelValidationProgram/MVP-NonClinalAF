# Description of output files

# File Descriptions

Below is a description of metadata for each set of output tables. The first number is the seed for the simulation.

Note: The origin of the simulations was to better understand how oysters adapt to salinity and temperature gradients in estuaries. Opt0 in the simulations corresponds to the salinity trait, which also corresponds to Env2 in the publications.


Shown is an example below for the seed `1231094`

## SLiM outputs _METADATA DONE, BUT NEEDS SECOND CHECK__
* `1231094.trees` The tree sequencing file that is input to pyslim. See the SLiM manual for details.
* `1231094_LA.txt` Local adaptation of the metapopulation through time. See Metadata for details.
* `1231094_fitnessmat.txt` an n_deme x n_deme table that indicates the mean fitness of individuals from the source deme (in columns) in the transplant deme (in rows)
* `1231094_fitnessmat_ind.txt` an n_deme x n_individual table that indicates the fitness of each source individual (in columns) in every transplant deme (in rows)
* `1231094_VCF_causal.vcf.gz` the VCF output from slim for the causal mutations - it contains some additional information not output by pyslim
	* Each row is a locus simulated in SLiM
	* Each column is an individual (all output)
	* The INFO field contains additional info about the mutation in SLiM that is not retained by pyslim
	* This table is redundant, but used to cross-check the pyslim outputs.
	* For more information on VCF formats, see http://samtools.github.io/hts-specs/
* `1231094_ind.txt` Individual information output last generation. See Metadata for details.
* `1231094_info.txt`  Simulation information. See Metadata for details.
* `1231094_mig_mat.txt` This is the migration matrix that was used in the sims. It is used in pyslim for recapitation. Without it, pyslim will run forever and never recapitate.
* `1231094_muts.txt` This is a table of information about each mutation that is simulated SLIM AND has an allele frequency > 0.01 or < 0.99. See Metadata for details.
* `1231094_outfile.error.txt` This is the command line error output `>>` for the SLiM run. 
* `1231094_outfile.txt` This is the command line standard output `>` for the SLiM run
* `1231094_popInfo.txt` deme information. See Metadata for details.
* `1231094_popInfo_m.txt` migration rates between populations, used for plotting. See Metadata for details.

## Output after pyslim and VCF filtering _METADATA DONE_
* `1231094_pytree.error.txt` This is the command line error output `>>` for the pyslim run
* `1231094_pytree.out.txt` This is the command line standard output `>>` for the pyslim run
* `1231094_plusneut_MAF01.recode2.vcf.gz`
	* The VCF file output after pyslim and filtered to an MAF of 0.01 across the entire population
	* Each row is a locus simulated in SLiM or added by pyslim
	* Each column is an individual (all output)
	* The ID of the causal mutations is in the ALT table and should match `mutID` in other tables
	* For more information on VCF formats, see http://samtools.github.io/hts-specs/


## Outputs after sampling and processing sims in R
* `1231094_R.error`  This is the command line error output `>>` for the R run
* `1231094_R.out` This is the command line standard output `>` for the R run
* `1231094_RDA.RDS` This is the R object output by the R function `rda()` run on the sampled individuals 
(10/deme x 100 demes = 1000 individuals), all mutations with MAF > 0.01, and the 2 simulated environments.
* `1231094_RDA_structcorr.RDS` This is the R object output by the R function `rda()`. This is a partial RDA run on the sampled individuals 
(10/deme x 100 demes = 1000 individuals), all mutations with MAF > 0.01, the 2 simulated environments, and including `Condition(PC1 + PC2)`, 
where PC1 and PC2 were the first two principal components of the population structure.
* `1231094_RDAcorrected_loadings.txt`
	* Loadings of each environmental variable on each RDA axis after structure correction, corresponding to `1231094_RDA_structcorr.RDS`
	* Row names correspond to `sal` (Env2) and `temp` selective environments
* `1231094_RDAloadings.txt`
	* Loadings of each environmental variable on each RDA axis, corresponding to `1231094_RDA.RDS`
	* Row names correspond to `sal` (Env2) and `temp` selective environments
* `1231094_Rout_Gmat_sample.txt` Genotypes from the sampled individuals (10/deme x 100 demes) at loci with MAF > 0.01 in the sample
	* Mutations are in rows and correspond to `1231094_Rout_muts_full.txt`
	* Individuals are in columns and correspond to `1231094_Rout_ind_subset.txt`
* `1231094_Rout_RDA_predictions.txt` correlation between (a linear prediction of the weighted RDA loadings) and (the individual's temp phenotype). See Metadata for more information.
* `1231094_Rout_af_pop.txt` Allele freq. of each population for each mutation after sampling of individuals (10 ind/deme x 100 demes)  and filtering MAF > 0.01. See Metadata for more information.
* `1231094_Rout_af_sal.txt` Allele freq. of all populations with the same salinity (Env2), for each mutation after sampling of individuals (10 ind/deme x 100 demes) and filtering MAF > 0.01. See Metadata for more information. 
* `1231094_Rout_af_temp.txt` Allele freq. of all populations with the same temperature, for each mutation after 
sampling of individuals (10 ind/deme x 100 demes) and filtering MAF > 0.01. See Metadata for more information.
* `1231094_Rout_ind_subset.txt` A table with information about each individual that was sampled for analysis (10 ind/deme x 100 demes)
See Metadata for more information.
* `1231094_Rout_muts_full.txt`  A table with information about each locus used for analysis (loci after sampling individuals
and then filtering MAF > 0.01). See Metadata for details
* `1231094_Rout_simSummary.txt`  See metadata in `summary_20220428.txt_metadata.md`
* `1231094_genotypes.lfmm` Transposed `1231094_Rout_Gmat_sample.txt` for running lfmm.
* `1231094_genotypes.pca` Written output of running `pca()` from `lfmm`. See https://academic.oup.com/mbe/article/36/4/852/5290100 for details.
* `1231094_genotypes.pcaProject` Written output of running `pca()` from `lfmm`. See https://academic.oup.com/mbe/article/36/4/852/5290100 for details.
* `1231094_lfmm2_sal.RDS` R object saved from the salinity (Env2) model with `lfmm2`
* `1231094_lfmm2_temp.RDS` R object saved from the temperature (Env2) model with `lfmm2`
* `1231094_pca.RDS` R object saved from the `pca()` function with `lfmm2`

## Figures related to processing sims in R
* `1231094_pdf_1pop.pdf` plots at the population level, fitness for each population
* `1231094_pdf_1pop_mig.pdf` migration patterns for the simulation
* `1231094_pdf_2muts.pdf` plots of individual mutation effect sizes
* `1231094_pdf_3bStructureCorrs.pdf` plots of how structure (PC axes) are related to environment or to LFMM structure estimates
* `1231094_pdf_3pca.pdf` principal components plots
* `1231094_pdf_4manhattan_LFMM.pdf` Manhattan plot for outliers from LFMM
* `1231094_pdf_5RDA.pdf` RDA plots for individuals or for mutations, with and without structure correction. Sometimes the same plot is twice but with different arrow scaling, I had trouble figuring out how to scale the arrows at first.
* `1231094_pdf_6manhattan_RDA.pdf` Manhattan plots for RDA outliers
* `1231094_pdf_6zRDA_predict.pdf` Accuracy of RDA prediction, with and without structure correction and for each trait
* `1231094_pdf_7_afcors.pdf` This plot shows the effect of pooling samples with the same environment together to calculate correlations between allele frequency and environment. Even very low correlations become artificially inflated.
* `1231094_pdf_8manhattan_cor.pdf` Manhattan plot for Kendall's rho cor(allele frequency, environment)
* `1231094_pdf_8manhattan_fst.pdf` Manhattan plot for FST calculated by OutFLANK (Weir and Cockerham)
* `1231094_pdf_8phen-env_af-env.pdf` visualization of phenotypic clines and allele frequency patterns for each simulated trait
* `1231094_pdf_heatmaps.pdf` genotype heatmaps for different subsamples of indivduals - subsampled by their environment or by their location on the landscape.
	* Individuals are in rows and causal loci are in columns. Yellow is 2 copies of the ancestral allele (0 effect on phenotype), Orange is heterozygote, and red is 2 copies of the derived allele (with a phenotypic effect)
