It's been a few months.

# First step was to sync git across my devices.

**On my computer:**

Ignored graphics files, committed changes to R code and synced to github 

**On Discovery:**

Ignored sim_output and graphics files, committed changes and synced to github

```
git status
git add #a shell file was changed and some conflicts in .gitignore were solved
git commit 
git pull # sync whatever is on github with discovery
git push
```

remember for a password I have to use my personal access token

# To Do List for sim setup

- On my computer, I recreated the list of simulation parameters, recall that these should produce architectures that don't take so long to coalesce
   - I also moved the files for the 20210920 runs to my archived folder

## To Do List for R code c-AnalyzeSimOutput.R

I initiated this with OOD, since all my simulation outputs are not synced to github. 

- [ ] add Va calc with and without filtered loci
   - [x] Based on the entire sample, currently I have:
      - `va_temp_full` and `va_sal_full` to describe the va of each mutation
      - `va_temp_total` and `va_sal_total` are the sums of the previous columns and also based on the entire sample
      - `va_temp_full_prop` and `va_sal_full_prop` the proportion of va based on the entire sample
      - [x] **DONE: add definitions to metadata**
   - [x] Based on the subsample, these columns have been called (not changing since it would affect too much code):
      - `Va_temp` and `Va_sal`
      - `Va_temp_prop` and `Va_sal_prop`
      - 
- Pop plots 
   - [x] the migration/Ne plot is too short - It was 5x4, I made it 5in x 5in
   - [x] remove numbers from the optimums plot
 - Phen-env Af-env plots
   -  [x] there seems to be a bug in the in the color ramp plotting legend (did I pick colors randomly?), 
      -    FIXED!!! this was around lines 1200 in the code.
   -  [x] bug the prop. VA plot fixed, and in the legend for prop. Va
 -  Number of neutral SNPs in demographies
   -    N equal\nm breaks producing twice as much as everything else
   -    N variable\nm variable produces on average the same number of SNPs, but sometimes produces 100,000 SNPs
   -    I have a vague memory of editing the sims to try to reduce the above effect, so for now nothing has been done. 

- [x] Calculate and plot FST for the sims 
   - [x]    muts_full$He_outflank   **ADDED TO METADATA**
   - [x]    muts_full$Fst_outflank   **ADDED TO METADATA**
   - [ ]    meanFst in sim output  **ADDED TO METADATA**
- [x] RDA
   - [x] `RDA1_propvar` and `RDA2_propvar` - **added to metadata**
   - [x] `subset_indPhen_df$RDA_predict_salPhen_20KSNPs` - this was previosuly called `t`
   - [x] FIX RDA MUTATION PREDICTION
      - [x] added `muts_full$RDA_mut_temp_pred` and `muts_full$RDA_mut_sal_pred`, need to **ADDED TO METADATA**
      - [x] I removed this code:

```  
  RDA_RDA1Loading_cor_tempEffect <- cor(muts_full$RDA1_score, muts_full$mutTempEffect, use="complete.obs")
  RDA_RDA1Loading_cor_salEffect <- cor(muts_full$RDA1_score, muts_full$mutSalEffect, use="complete.obs")
  RDA_RDA2Loading_cor_tempEffect <- cor(muts_full$RDA2_score, muts_full$mutTempEffect, use="complete.obs")
  RDA_RDA2Loading_cor_salEffect <- cor(muts_full$RDA2_score, muts_full$mutSalEffect, use="complete.obs") 
```
      - [x] and replaced it with this code: 
```
RDA_RDALoading_cor_tempEffect <- cor(muts_full$RDA_mut_temp_pred, muts_full$mutTempEffect, use="complete.obs")  **ADDED TO METADATA**
RDA_RDALoading_cor_salEffect <- cor(muts_full$RDA_mut_sal_pred, muts_full$mutSalEffect, use="complete.obs")  **ADDED TO METADATA**
```

   - [x] output: "X....RDA1_temp_cor" - I think I fixed this
   - [x] RDA plot too busy - black diamonds


# Setting up for 20220117 Run
* Starting with the first 1000 fast simulation parameters on lotterhos partiition. If I understand the architecture correctly, each node has 36 cores and each core has 5 GB. So I should be able to run: 5 GB/core * 36 cores * 2 nodes / 2.5GB per job = 144 jobs at a time. I'm going to try:
```
#SBATCH --partition=lotterhos
#SBATCH --mem=2G
#SBATCH --nodes=1
#SBATCH --array=2-1000%70
```

`sbatch d-run_nonAF_sims_0Slim-fastruns-20220117.sh`

Submitted batch job 22811525

`squeue -u lotterhos` to check the jobs of the user

`squeue -p lotterhos` to check what's running on our nodes

Edit the R script to update with the date of the runs `c-AnalyzeSimOutput.R` on this line:
`allsims <- load("src/0b-final_params-20220117.RData")`

Run the R code for the first 1000 runs:

`sbatch src/e-run_nonAF_sims_1R-fastruns-20220117.sh`

Submitted batch job 22834818

## First 1000 runs - did they work?
```
ls -l *_muts.txt | wc -l # this is the number of sims that completed the SliM run
999

ls -l *_plusneut_MAF01.recode2.vcf.gz | wc -l # this is the number of sims that completed the vcf filtering
999

ls -l *_pdf_1pop.pdf | wc -l # number of sims that were analyzed through the population step
984

ls -l *_2muts.pdf | wc -l # number that analyzed the mutations
984

ls -l ../*.pcaProject | wc -l # this is the number of sims that were analyzed through the PCA step in the R script
804

ls -l *_Rout_simSummary.txt | wc -l # this is the number of sims that were analyzed through the final output in the R script
701
```

### That's weird

FIGURE OUT WHAT THE HELL IS GOING ON! From the numbers above, I can tell something is happening in the R script. It's not finishing on all simulations.

## In the mean time, I think I can run the rest of the simulations.

`sbatch d-run_nonAF_sims_0Slim-fastruns-20220117-b.sh` after updating the file in the shell script to run the correct params

Submitted batch job 22887383

```
ls -l *_muts.txt | wc -l # this is the number of sims that completed the SliM run
1950 # all ran!

ls -l *_plusneut_MAF01.recode2.vcf.gz | wc -l 
1950
```

## R code to do

-[ ] Figure out which simulations aren't finishing the beginning of the R script and which aren't getting to the end of the R script

OK, so one sim that didn't finish was 1231097. The sim failed at the Outflank stage, because some markers had NAs. In tracing the issue, I 
found that some loci had triallelic sites (sites should be a combo of 0 and 1):

```
table(vcf_full@gt[8993,])

 0|0  0|1  0|2  1|0  1|1  2|0  2|2   GT 
8802  173  318  202   60  384   61    1 
```

After filtering, this was a single site in the genome:

```
which(!complete.cases(G_full_subset))
8869
```

Triallelic sites are typically removed from real data, and they are rare in the simulations, so removing these sites should not affect the results.
To track this, I stored the number of triallelic sites in the variable `num_multiallelic`

## Rerunning R code

```
sbatch src/e-run_nonAF_sims_1R-fastruns-20220117.sh
Submitted batch job 22914512
```

## Graphs TO DO:
Previous graphs can be found at `/work/lotterhos/MVP-NonClinalAF/run20210920/20210920_graphs/`

-  [ ] Put together a 10x2 figure for the SS and Est demographies
-  [ ] Fst across sims
-  [ ] The phenotypic clines, AF clines, and heatmap figure
-  [ ] LFMM figures - what to show?
-  [ ] RDA figure 1 - demography histograms, demography scatterplots
-  [ ] RDA figure 2 - a 6 panel plot - SS Clines, SS Mtn, and Est - mutation loading on top, individual scores on bottom
- [ ] create main text figure for simple scenarios
- [ ] create a supplemental figure for the complex scenarios

Potential additional steps:
- [ ] run outflank for FST outliers
- [ ] run higher gene flow sims with less LA in SS scenarios
- [ ] add fake correlated environments to RDA analysis
- [ ] add ROC scores to output
- [ ] add sims that have demog correlated with a non-selective environment

## Housekeeping
- [ ] Download YML files from /work/lotterhos/MVP-NonClinalAF/src to https://github.com/northeastern-rc/lotterhos_group
