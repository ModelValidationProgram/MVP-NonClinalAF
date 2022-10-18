## `summary_LAthroughtime_20220428_20220726.txt` 

A table summarizing the local adaptation through time in each simulation.

This data was output from the SLiM run and compiled for analysis. The data was used
to check that the amount of local adaptation reached equilibrium.

-------------------------
### Basic outputs
-------------------------
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