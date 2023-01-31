`/multipheno_multienv/src` folder

- `a-BIOCLIM_extraction.R` This script extracts BIOCLIM variables for any given lat/long grid and prepares files for input to SLiM
- `b-non_wf_range_exp_working.slim` the SLIM code for multivariate range expansion with 6 traits in 6 environments in continuous space
- `c-pyslim.py` pyslim recapitation script for the `b-non_wf_range_exp_working.slim`
- `c-pyslim.sh` shell script to run `c-pyslim.py` on NU Discovery Cluster
- `d-processVCF.Rmd` R code for sampling 1000 individuals from the SLiM output and running the various RDA models