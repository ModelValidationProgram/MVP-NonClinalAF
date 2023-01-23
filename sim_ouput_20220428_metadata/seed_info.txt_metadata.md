Most of this information is redundant with the parameters put into the simulation.
These parameters are documented in the `/src/0a-setUpSims.Rmd` code and `0b-final_params-20220428.txt`

### seed_info.txt (SLiM output) Simulation information

* `seed`  simulation seed
* `MIG_x` migration across longitude populations
* `MIG_y` migration across latitude populations
* `mu` mutation rate
* `r` recombination rate
* `BURNIN_1` length of first burnin (typically 0 opimums across metapopulation)
* `BURNIN_2` length of second burnin (typically a transition from 0 opimums across metapopulation to final environment values)
* `TOTAL_GEN` the total number of generations
* `METAPOP_SIDE_x` the number of metapopulations across longitude
* `METAPOP_SIDE_y` the number of metapopulations across latitude
* `SIGMA_K_0` Strength of stabilizing selection on first trait
* `SIGMA_K_1` Strength of stabilizing selection on second trait
* `SIGMA_K_Cov` Covariance in selection on both traits (always 0 in this study)
* `SIGMA_QTN_1` Standard deviation of normal distribution from which effect sizes are drawn for first trait
* `SIGMA_QTN_2` Standard deviation of normal distribution from which effect sizes are drawn for second trait
* `SIGMA_QTN_Cov` Covariance in multivariate-normal distribution from which effect sizes are drawn for both traits (always 0 in this study)
* `demog` Landscape type for simulation
* `Ntraits` Number of traits - 1 or 2
* `ispleiotropy` 0 for no pleiotropy, 1 for pleiotropy
* `iscontrol` If true, all individuals have equal fitness. (always false in this study)
* `Nequal` Setting corresponds to different patterns of N on the landscape. See main paper for details.
* `isVariableM` If true, migration rates between populations are chosen from a distribution.
* `MIG_breaks` If true, migration breaks occur at two latitudes. See main paper for details.