`multipheno_multienv/output_multisim` folder

Outputs for the 6-trait 6-environment continuous space multivariate simulations

- `_pytree.error.txt` output file for errors from running pyslim (empty)
- `_pytree.out.txt` output from pyslim

- `892657863_ind.txt` Information about each simulated individual at the end of the sim
	- `seed` simulation seed
	- `ind_index` index of of the individual
	- `x`  x location on the landscape
	- `y`  y location on the landscape
	- `phenotype1_mat`  individual trait value for this environmental trait
	- `phenotype2_MTWetQ`  individual trait value for this environmental trait
	- `phenotype3_MTDQ`  individual trait value for this environmental trait
	- `phenotype4_PDM`  individual trait value for this environmental trait
	- `phenotype5_PwarmQ`  individual trait value for this environmental trait
	- `phenotype6_PWM`  individual trait value for this environmental trait
	- `env1_mat`  environment where the individual was sampled
	- `env2_MTWetQ`  environment where the individual was sampled
	- `env3_MTDQ`  environment where the individual was sampled
	- `env4_PDM`  environment where the individual was sampled
	- `env5_PwarmQ`  environment where the individual was sampled
	- `env6_PWM` environment where the individual was sampled

- `892657863_muts.txt` Information about each QTN at the end of the sim
	- `seed`  simulation seed
	- `mutID`  ID of the mutation
	- `muttype`  mutation type in SliM
	- `p`  allele frequency at the end of the simulation across all individuals
	- `mut1-mat-Effect`  effect size of derived allele on this environmental trait
	- `mut2-MTWetQ-Effect`  effect size of derived allele on this environmental trait
	- `mut3_MTDQ-Effect`  effect size of derived allele on this environmental trait
	- `mut4_PDM-Effect`  effect size of derived allele on this environmental trait
	- `mut5_PwarmQ-Effect`  effect size of derived allele on this environmental trait
	- `mut6_PWM-Effect`  effect size of derived allele on this environmental trait

- `892657863_pca.RDS` an R data object for the output of the principal components analysis

- `892657863_plusneut_MAF01.recode2.vcf.gz` the vcf file for all individuals including neutral mutations output after tree sequencing and filtering for MAF > 0.01

- `892657863_simulation_parameters.txt` parameters that correspond to the run of `../src/b-non_wf_range_exp_working.slim`


- `892657863_throughtime.txt` Information about the simulation through time
	- `sim.generation`	 generation
	- `num_muts`	number of QTN mutations
	- `num_ind`	  number of individuals
	- `mean-phenotype1-mat`	 mean trait value across all individuals for this environmental trait
	- `mean-phenotype2-MTWetQ`	 	 mean trait value across all individuals for this environmental trait
	- `mean-phenotype3-MTDQ`		 mean trait value across all individuals for this environmental trait
	- `mean-phenotype4-PDM`	 	 mean trait value across all individuals for this environmental trait
	- `mean-phenotype5-PwarmQ`		 mean trait value across all individuals for this environmental trait 
	- `mean-phenotype6-PWM`	 	 mean trait value across all individuals for this environmental trait
	- `sd-phenotype1-mat`	 standard deviation in trait value across all individuals for this environmental trait
	- `sd-phenotype2-MTWetQ`standard deviation in trait value across all individuals for this environmental trait	 
	- `sd-phenotype3-MTDQ`	standard deviation in trait value across all individuals for this environmental trait 
	- `sd-phenotype4-PDM`	standard deviation in trait value across all individuals for this environmental trait 
	- `sd-phenotype5-PwarmQ` standard deviation in trait value across all individuals for this environmental trait	 
	- `sd-phenotype6-PWM`	standard deviation in trait value across all individuals for this environmental trait 
	- `corr_phen_env_1mat`	 correlation between this environmental trait and the environmental variable
	- `corr_phen_env_2MTWetQ`	correlation between this environmental trait and the environmental variable 
	- `corr_phen_env_3MTDQ`	 correlation between this environmental trait and the environmental variable
	- `corr_phen_env_4PDM`	 correlation between this environmental trait and the environmental variable
	- `corr_phen_env_5PwarmQ`	correlation between this environmental trait and the environmental variable 
	- `corr_phen_env_6PWM`  correlation between this environmental trait and the environmental variable

- `892657863_VCF_causal.vcf`
	- VCF file output by SliM at the end of the simulation. The file contains information about each QTN that is lost after recapitation in pyslim.


- `892657863.trees` trees file output by pyslim


- `892657863genotypes.txt` Genotype matrix for the 1000 sampled individuals in rows and SNPs in columns.
