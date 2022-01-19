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

Submitted batch job 22823441

## Graphs TO DO:
Previous graphs can be found at `/work/lotterhos/MVP-NonClinalAF/run20210920/20210920_graphs/`

-  [ ] Put together a 10x2 figure for the SS and Est demographies
-  [ ] Fst across sims
-  [ ] The phenotypic clines, AF clines, and heatmap figure
-  [ ] LFMM figures - what to show?
-  [ ] RDA figure 1 - demography histograms, demography scatterplots
-  [ ] RDA figure 2 - a 6 panel plot - SS Clines, SS Mtn, and Est - mutation loading on top, individual scores on bottom
  - create main text figure for simple scenarios
  - create a supplemental figure for the complex scenarios

Potential additional steps:
- [ ] run outflank for FST outliers
- [ ] run higher gene flow sims with less LA in SS scenarios
- [ ] add fake correlated environments to RDA analysis
- [ ] add ROC scores to output

## Housekeeping
- [ ] Download YML files from /work/lotterhos/MVP-NonClinalAF/src to https://github.com/northeastern-rc/lotterhos_group
