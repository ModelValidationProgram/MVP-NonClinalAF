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
-  `1231094_plusneut_MAF01.recode2.vcf.gz`
-  * `1231094_Rout_Gmat_sample.txt` Genotypes from the sampled individuals (10/deme x 100 demes) at loci with MAF > 0.01 in the sample
	* Mutations are in rows and correspond to `1231094_Rout_muts_full.txt`
	* Individuals are in columns and correspond to `1231094_Rout_ind_subset.txt`
- * `1231094_plusneut_MAF01.recode2.vcf.gz`
  - * `1231094_popInfo.txt` deme information. See Metadata for details.
- * `1231094_muts.txt` This is a table of information about each mutation that is simulated SLIM AND has an allele frequency > 0.01 or < 0.99. See Metadata for details.
- * `1231094_VCF_causal.vcf.gz` the VCF output from slim for the causal mutations - it contains some additional information not output by pyslim
	* Each row is a locus simulated in SLiM
	* Each column is an individual (all output)
	* The INFO field contains additional info about the mutation in SLiM that is not retained by pyslim
	* This table is redundant, but used to cross-check the pyslim outputs.
	* For more information on VCF formats, see http://samtools.github.io/hts-specs/

* `1231094_ind.txt` Individual information output last generation. See Metadata for details.
* `1231094_info.txt`  Simulation information. See Metadata for details.

- ALL PDF FILES
- 
