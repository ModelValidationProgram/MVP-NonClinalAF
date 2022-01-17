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

## To Do List for R code
- [ ] add Va calc with and without filtered loci
- Pop plots 
   - [ ] the migration/Ne plot is too short
   - [ ] remove numbers from the optimums plot
 - Phen-env Af-env plots
   -  [ ] there seems to be a bug in the in the color ramp plotting legend (did I pick colors randomly?), the prop. VA plot, and in the legend for prop. Va
 -  Number of neutral SNPs in demographies
   -    N equal\nm breaks producing twice as much as everything else
   -    N variable\nm variable produces on average the same number of SNPs, but sometimes produces 100,000 SNPs
- [ ] output: "X....RDA1_temp_cor"
- [ ] RDA plot too busy - black diamonds
- [ ] Calculate FST for the sims 
- [ ] FIX RDA MUTATION PREDICTION

# Graphs TO DO:
-  [ ] Put together a 10x2 figure for the SS and Est demographies
-  [ ] The phenotypic clines, AF clines, and heatmap figure
-  [ ] LFMM figures - what to show?
-  [ ] RDA figure 1 - demography histograms, demography scatterplots
-  [ ] RDA figure 2 - a 6 panel plot - SS Clines, SS Mtn, and Est - mutation loading on top, individual scores on bottom
  - create main text figure for simple scenarios
  - create a supplemental figure for the complex scenarios
