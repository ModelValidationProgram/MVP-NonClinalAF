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

## Sync files to git
```
git status
git add <added files in commit stages>
git commit
git pull
git push
```

## Rerunning R code

Before rerunning, I deleted some of the R outputs, so I can be sure the files are replaced:

```
cd sim_output_20220428/
rm *.pdf
rm *_Rout_simSummary.txt
rm -r  *genotypes.pca
```

I will also check the dates on the replaced files to make sure they are for the most recent run.

```
sbatch e-run_nonAF_sims_1R-fastruns-20220428.sh #Submitted batch job 32304308

sbatch e-run_nonAF_sims_1R-fastruns-20220428-b.sh #Submitted batch job 32500214

sbatch e-run_nonAF_sims_1R-longruns-20220428.sh #Submitted batch job 32320018
```

Check runs are finished (first and last outputs of R script)
```
ls -la *_pdf_1pop.pdf | grep "Oct" | wc -l #should be 2250

ls -la *_Rout_simSummary.txt | grep "Oct" | wc -l #should be 2250
```

additional checks (I overwrote these files, so just want to check all the files are recent and there aren't any old ones left)

```
cd sim_output_20220428
ls -la | head -n 100 # check dates on files
ls -la | grep "1231094" | wc -l #number of files expected for each sim = 50 files
ls -la | grep "July" | wc -l # should be ? based on Slim runs
ls -la | grep "Oct" | wc -l # should be ? based on rerun of R code

```
## Wed Nov 9 Check runs are complete

ls sim_output_20220428/*_2muts.pdf | wc -l

ls -l *_Rout_simSummary.txt | wc -l 

2250 Simulations are complete!

## Sat Nov 12
Cat outputs together. I accidently overwrote the 20220726 outputs - so those are the final files


## Outputs to copy for Jeff Ross-Ibarra

-[x] https://github.com/ModelValidationProgram/MVP-NonClinalAF/blob/main/notebook/20220914_JeffRossIbarra.md

## Multvariate Sims
I dug up an old 2-trait range expansion continuous space simulation from TTT. Edits:
- used 6 standardized WorldClim variables that were weakly correlated with each other
- decided on a system to include all 6 traits in each simulation; whether mutations effect a trait depends on the parameters. 
  - I didn't notice a major performance difference using this system, as the unselected traits are just 0
-  I tested the sims from 2 traits to 6 traits and they all seem to be working well!

**Oct 24 2022**
- Running recapitation and vcf filtering on multivariate simulations - job ID 32337815

- Running recapitation and vcf filtering on multivariate simulations - job ID 32337835

**Oct 27 2022**
- Bug in getting the individual environments in the output. It shouldn't have affected the trees file, but running everything again just in case.

