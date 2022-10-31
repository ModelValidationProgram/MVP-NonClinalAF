
## Using old TTT code
`TTT_2pheno_2envi` renamed ---> `MVP_multipheno_multienvi`

### To do in slim / on computer:
- [x] check SLiM version is the same as OOD
- [x] update to dmv_norm. I tried this and the behavior was unexpected, it works without it so will skip for now.
- [x] output mutation effect size table
- [x] output individual environments and phenotypes table
- [x] output mutation vcf file
- [x] compare parameters to the MVP simulation
- [x] use tree sequencing to recap and output VCF - see Recombination paper

### In R:
- [x] run RDA and test ideas with the 2 trait sims
- [x] evaluate clinal paradigm with proportion of causal clines
- [x] add non-adaptive environments to RDA and test it

## Expand sim to multiple traits
- [x] for simplicity just have one mutation type and pleiotropic mutations, no trade-offs


## See timestamps in last notebook post
After fixing the errors, I ran it again just in case:
```
(base) [lotterhos@ood src]$ sbatch c-pyslim.sh
Submitted batch job 32368804
```

## Base simulations
- [ ] Remake manuscript figures
- [ ] Remake/Rearrange Figure 5
- [ ] Conduct GWAS analysis
- [ ] outputs for Jeff R-I

## Multitrait simulations
- [x] Run everything from the start in R, to make sure it works
- [x] clean up code
- [x] Run a PCA just to visualize structure, colored by xy location

- [x] Main plot - 6 environment panels, PCA, and a barplot showing accuracy of the RDA prediction in the 3 scenarios (or 4, if add PCA)
- [x] Supp plot - correlation matrices for bioclim environments
- [x] Supp plot - correlation matrices for mutation effect size
- [x] Supp plot - RDA colored by each environmental variable (for the base RDA)
- [x] write subsetted datasets to file
- [ ] create a tutorial
- [ ] edit movie to make 1 min long (7x speed)
