Working on responding to reviewer comments

- [x] Create a new figure from SS Mtn to show that different architectures can plot to the same phenotype in RDA space. `g4-DifferentArchSameGeno.Rmd`
- [x] standardized phenotype and mutation z-score (use `scale()` function)
- [ ] Framework for determining clinal paradigm - Trait level - using GWAS and seeing if GWAS hits show clines or not

## Framework for assessing clinal paradigm

New outputs:
These two are the true positive rate of GWAS hits (proportion of Benjamini-Hotchberg corrected P-values < 0.05):
- `gwas_TPR_sal`
- `gwas_TPR_temp`

The next two are false discovery rate of GWAS hits (proportion of Benjamini-Hotchberg corrected P-values < 0.05), excluding neutral-linked loci affected by selection:
- `gwas_FDR_sal_neutbase`
- `gwas_FDR_temp_neutbase`  

These two are the proportion of the top 5% of GWAS hits that show clines in allele frequency without structure correction:
- `clinalparadigm_sal_proptop5GWASclines`
- `clinalparadigm_temp_proptop5GWASclines`
This works for the simulation I tried it on in the moderately polygenic case. (I tried using LFMM outliers to estimate clines, but it didn't work for temperature because of overcorrection for structure)

These two are the proportion of the significant GWAS hits that show clines in allele frequency without structure correction:
- `clinalparadigm_sal_propsigGWASclines`
- `clinalparadigm_temp_propsigGWASclines`
