This file provides all the information that is used to setup the simulations in SLiM.

For details see `0a-setUpSims.Rmd`, which is used to create this file.

* `level` A full description of the sim, including the architecture, landscape, and demography simulated. 225 levels
* `reps` The replicate number. Each level has 10 replicates
* `arch` A description of the genetic architecture (15 levels)
* `demog_name` A description of the landscape-demography (15 levels)
* `demog_level_sub` A description of the demography
* `demog_level` A description of the landscape
* `MIG_x` Migration rate along the x-axis (unless otherwise specified)
* `MIG_y` Migration rate along the y-axis (unless otherwise specified by the landscape or `isVariableM`)
* `xcline` describes the relationship between the x location and environment - `linear` or `V`  
* `ycline` describes the relationship between the y location and environment - always `linear` 
* `demog` describes the landscape - `SS` for stepping stone or `Estuary`
* `METAPOP_SIDE_x` number of demes on the x axis - always 10
* `METAPOP_SIDE_y` number of demes on the y axis - always 10
* `Nequal` an indicator for the carrying capacity ($N_deme$) for each deme on the landscape. This is implemented with custom code in SLiM.
	* `0` equal $N_deme$ across the landscape
	* `2` a North-to-South cline in $N_deme$
	* `3` $N_deme$ highest in the center and lowest at the edges
	* `4` variable $N_deme$
* `isVariableM` an indicator as to whether migration rate is constant or variable. This is implemented with custom code in SLiM.
	* `0` use the migration rate specified by MIG_x and MIG_y
	* `1` use variable migration rates implemented with custom code in SLiM
* `MIG_breaks` an indicator as to whether migration breaks are added on the landscape. This is implemented with custom code in SLiM.
	* `0` do not implement custom code
	* `1` implement custom code
* `arch_level_sub` The five genetic architecture levels simulated within each `arch_level`
 	* "1-trait", "2-trait-no-pleiotropy-equal-S", "2-trait-no-pleiotropy-unequal-S", "2-trait-pleiotropy-equal-S", "2-trait-pleiotropy-unequal-S"
* `arch_level` The three genic levels: "oligogenic", "mod-polygenic", "highly-polygenic"
* `MU_base` The baseline mutation rate for neutral loci 1e-07 
* `MU_QTL_proportion` The mutation rate for QTN loci = `MU_base` x `MU_QTL_proportion`
* `SIGMA_QTN_1`  Standard deviation in distribution of effect size of a new mutation on trait 1 (x axis trait - Env2)
* `SIGMA_QTN_2` Standard deviation in distribution of effect size of a new mutation on trait 2 (y axis trait - temp)
* `SIGMA_K_1` Standard deviation in fitness function - strength of stabilizing selection - on trait 1 (x axis trait - Env2)
* `SIGMA_K_2` Standard deviation in fitness function - strength of stabilizing selection - on trait 2 (y axis trait - temp)
* `N_traits` Number of traits simulated (1 for temp trait only; 2 for temp and Env2 trait)
* `ispleiotropy` Whether or not plietropy was modeled. 
	* If `0`, every mutation could only effect 1 trait.  Mutation effect size drawn from a normal distribution for the trait.
	* If `2`, every mutation could effect both traits. Mutation effect size drawn from a multivariate normal distribution for the traits.
* `seed` simulation ID and also a seed for SLiM
