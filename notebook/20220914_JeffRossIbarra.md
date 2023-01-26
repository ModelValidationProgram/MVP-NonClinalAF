# People

Jeff Ross-Ibarra
Brandon 
Katie

# Agenda

- Plan for future use of sims: Lotterhos Lab
- Plan for future use of sims: Ross-Ibarra Lab
  - Breeding simulations 
- Structure of simulation code
- Plan for achival
  - Plan for archival is with BCO-DMO after the reviews come back and a 1 year embargo
  - This plan is in place to give Brandon time to publish additional papers with the simulations
- Discussion of outputs Ross-Ibarra might want


Outputs:

PRE-SAMPLING (all ~10,000 individuals)
-  `1231094_plusneut_MAF01.recode2.vcf.gz` vcf file of genotypes for all ~10,000 individuals in the simulation
	- `1231094_ind.txt` Individual information output last generation.
	- `1231094_muts.txt` This is a table of information about each QTN that is simulated SLIM AND has an allele frequency > 0.01 or < 0.99. 
- `1231094_popInfo.txt` deme information. 
- `1231094_VCF_causal.vcf.gz` the VCF output from slim for the causal mutations - it contains some additional information not output by pyslim, as well as rare alleles that are filtered with the MAF
	* Each row is a locus simulated in SLiM
	* Each column is an individual (all output)
	* The INFO field contains additional info about the mutation in SLiM that is not retained by pyslim
	* This table is redundant, but used to cross-check the pyslim outputs.
- `1231094_info.txt`  Parameters for this seed


POST-SAMPLING AND ANALYSIS (sample of 1000 individuals from the landscape)
-  `1231094_Rout_Gmat_sample.txt` Genotypes from the sampled individuals (10/deme x 100 demes) at loci with MAF > 0.01 in the sample
	* Mutations are in rows and correspond to `1231094_Rout_muts_full.txt`
	* Individuals are in columns and correspond to `1231094_Rout_ind_subset.txt`
- ALL PDF FILES visualizing outliers and population structure
