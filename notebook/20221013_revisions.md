Working on responding to reviewer comments

- [x] Create a new figure from SS Mtn to show that different architectures can plot to the same phenotype in RDA space. `g4-DifferentArchSameGeno.Rmd`
- [x] standardized phenotype and mutation z-score (use `scale()` function)
- [ ] Framework for determining clinal paradigm - Trait level - using GWAS and seeing if GWAS hits show clines or not

## Framework for assessing clinal paradigm

New outputs:
These two are the true positive rate of GWAS hits (proportion of Benjamini-Hotchberg corrected P-values < 0.05):
- `gwas_TPR_sal`
- `gwas_TPR_temp`

These two are the proportion of the top 5% of GWAS hits that show clines in allele frequency without structure correction:
- `clinalparadigm_sal_propGWASclines`
- `clinalparadigm_temp_propGWASclines`
This works for the simulation I tried it on in the moderately polygenic case. (I tried using LFMM outliers to estimate clines, but it didn't work for temperature because of overcorrection for structure)

