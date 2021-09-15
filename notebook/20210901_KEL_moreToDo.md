# Aug 2021

### To Do: figure out why these sims take so long: `Est-Clines_N-equal_m_breaks` `Est-Clines_N-variable_m-variable`

One of them finished! But it took a day. Others are still running.

*oliogenic_1-trait__Est-Clines_N-equal_m_breaks*

```
(MVP_env) [lotterhos@login-01 src]$ seff 20883113_4
Job ID: 20883116
Array Job ID: 20883113_4
Cluster: discovery
User/Group: lotterhos/users
State: COMPLETED (exit code 0)
Cores: 1
CPU Utilized: 1-00:05:37
CPU Efficiency: 99.62% of 1-00:11:06 core-walltime
Job Wall-clock time: 1-00:11:06
Memory Utilized: 1.33 GB
Memory Efficiency: 66.69% of 2.00 GB
```

*polygenic_2-trait-no-pleiotropy-unequal-S__Est-Clines_N-variable_m-variable*

```
(MVP_env) [lotterhos@login-01 src]$ seff 20883113_111
Job ID: 20883904
Array Job ID: 20883113_111
Cluster: discovery
User/Group: lotterhos/users
State: COMPLETED (exit code 0)
Cores: 1
CPU Utilized: 5-04:36:12
CPU Efficiency: 98.23% of 5-06:51:09 core-walltime
Job Wall-clock time: 5-06:51:09
Memory Utilized: 18.76 GB
Memory Efficiency: 938.05% of 2.00 GB
```

*oliogenic_2-trait-no-pleiotropy-equal-S__Est-Clines_N-equal_m_breaks*
```
Job ID: 20883917
Array Job ID: 20883113_119
Cluster: discovery
User/Group: lotterhos/users
State: COMPLETED (exit code 0)
Cores: 1
CPU Utilized: 01:02:39
CPU Efficiency: 99.18% of 01:03:10 core-walltime
Job Wall-clock time: 01:03:10
Memory Utilized: 1.60 GB
Memory Efficiency: 79.88% of 2.00 GB
```

- The rest of them are still running after 6 days:
```
- (MVP_env) [lotterhos@login-01 src]$ squeue -u lotterhos
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
      20883113_141 lotterhos SlimRun2 lotterho  R 6-12:07:45      1 d3037
       20883113_96 lotterhos SlimRun2 lotterho  R 6-13:53:52      1 d3037
       20883113_81 lotterhos SlimRun2 lotterho  R 6-14:14:00      1 d3037
       20883113_66 lotterhos SlimRun2 lotterho  R 6-14:42:44      1 d3037
       20883113_51 lotterhos SlimRun2 lotterho  R 6-15:10:44      1 d3037
          20919505 lotterhos     bash lotterho  R 4-16:53:16      1 d3037
        20883113_6 lotterhos SlimRun2 lotterho  R 6-15:41:45      1 d3037
       20883113_36 lotterhos SlimRun2 lotterho  R 6-15:41:45      1 d3037
```

*20883113_141 polygenic_2-trait-pleiotropy-unequal-S__Est-Clines_N-variable_m-variable*


This simulation was cut off at 7 days. It produced a trees file with 118 million mutations (yikes). The "Plus Neuts" file was 3.92 tb, and after MAF filtering it was 2.39 TB. Yikes!
```
(base) [lotterhos@login-00 ~]$ seff 20883113_141
Job ID: 20883977
Array Job ID: 20883113_141
Cluster: discovery
User/Group: lotterhos/users
State: TIMEOUT (exit code 0)
Cores: 1
CPU Utilized: 6-22:16:18
CPU Efficiency: 98.97% of 7-00:00:28 core-walltime
Job Wall-clock time: 7-00:00:28
Memory Utilized: 72.60 GB
Memory Efficiency: 3630.16% of 2.00 GB
```


*20883113_96 polygenic_2-trait-no-pleiotropy-equal-S__Est-Clines_N-variable_m-variable*

This is seed 1231188. This was also cut off after 7 days. It doesn't look like it was finished recapitating.

```
Job ID: 20883832
Array Job ID: 20883113_96
Cluster: discovery
User/Group: lotterhos/users
State: TIMEOUT (exit code 0)
Cores: 1
CPU Utilized: 7-00:00:06
CPU Efficiency: 100.00% of 7-00:00:29 core-walltime
Job Wall-clock time: 7-00:00:29
Memory Utilized: 1.14 GB
Memory Efficiency: 57.02% of 2.00 GB
```


#### Create new slim file a-PleiotropyDemog_202109.slim

The above issue seems to be mostly with the Est-Clines_N-variable_m-variable demography, which needed a very large amount of memory and took a long time to run.

Migration rates were sampled from `sample(c(rep(0.001,6),0.001,0.01,0.1,0.25), 1)`, which I changed to `sample(c(0.001,0.01,0.1,0.25))`

I kept the sampling for the SS to `sample(c(rep(0.001,6),0.01,0.1,0.25)`, since this seemed to create cool patterns


#### Trial run d-run_nonAF_sims_0Slim_20210909.sh with  a-PleiotropyDemog_202109.slim and parameters that were tricky before

`sbatch d-run_nonAF_sims_0Slim_20210909.sh` - set it to line run 141, which was the one above that produced really large files.

Submitted batch job 21016929

The VCF file produced is still huge, > 1 TB, it takes a long time to filter the MAF. done in 4 days!

`seff 21016929`:

```
Job ID: 21016929
Array Job ID: 21016929_141
Cluster: discovery
User/Group: lotterhos/users
State: COMPLETED (exit code 0)
Cores: 1
CPU Utilized: 3-10:02:19
CPU Efficiency: 98.30% of 3-11:27:27 core-walltime
Job Wall-clock time: 3-11:27:27
Memory Utilized: 33.74 GB
Memory Efficiency: 337.40% of 10.00 GB
```

Ran in R:
(base) [lotterhos@login-01 src]$ sbatch e-run_nonAF_sims_1R_20210909.sh
Submitted batch job 21062450

taking hours to read in VCF file!

## To Do

- [x] While I'm waiting for the slow sims to finish:
  - [x] **metadata for output dataframe** https://github.com/ModelValidationProgram/MVP-NonClinalAF/blob/alan/sim_outputs_meta/README.md
  - [x] sync github
  - [x] check outputs for other parameter levels and write a notebook post (e.g. m breaks scenarios)
  - [x] Work on methods for manuscript
  - [x] **Develop full analysis script**
    - [x] List of questions
    - [x] Brainstorm figures to answer the questions
    - [x] cat together all the outputs
    - [x] write R code to analyze and make figures


- [ ] **R code** 
  - [x] Errors
    - [x] `Error in hist.default(muts_full$mutTempEffect, xlim = c(-1, 1), xlab = "Mutation temp effect",  : 
  some 'x' not counted; maybe 'breaks' do not span range of 'x'
Calls: hist -> hist.default`
    - [x] `Warning message:
position_dodge requires non-overlapping x intervals 
Warning message:
position_dodge requires non-overlapping x intervals 
Error in mod_temp@U[, 2] : subscript out of bounds`

  - [x] set salinity bubbles to be dark grey, decrease size of graphic output
  - [x] export RDA pc loadings for individuals. I bet this is interesting for the m-variable case.
  - [x] should I include the low MAF causal loci in the analysis? what about the highly polygenic case?
      - [x] I took these out. I was too worried about how the polygenic ones might bias results. Number removed is in `numCausalLowMAFsample`  
      - [x] `num_causal_postfilter` is the number of causal after sampling and MAF filtering of sample
      - [x] `num_causa_prefilter` is the number of causal loci in SliM at end of sim

  - [x] plot VA-prop vs. af cline
  - [x] change background for mutation AF clines density plot and explore different histogram types (density - what I have- is the best)
  - [x] COR (MUT EFFECT ON PHENOTYPE, MUT LOADING ON RDA AXIS) - test on an unequal S scenario
  - [x] remove row names from output table
  - [x] output number of outliers for each analysis


geom_point(aes(x=x, y=y,size=opt0), color="grey30") + geom_point(aes(x=x, y=y, color=opt1), size=2.5)

I change the RDA Qvalue to be consistent with lfmm cutoff

Variables added to analysis:
* `cor_af_temp_noutliers` number of outliers for cor(af,temp) after Bonferroni correction
* `cor_af_sal_noutliers` number of outliers for cor(af,salinity) after Bonferroni correction
* `nSNPs` total number of SNPs in analysis
* `cor_FPR_temp_neutSNPs` false positive rate in cor(af,temp) after Bonferroni correction, based on loci unaffected by selection
* `cor_FPR_sal_neutSNPs` false positive rate in cor(af,sal) after Bonferroni correction, based on loci unaffected by selection
* `LEA3.2_lfmm2_mlog10P_tempenv_noutliers` number of outliers for the lfmm temp model (qvalue <0.05)
* `LEA3.2_lfmm2_mlog10P_salenv_noutliers` number of outliers for the lfmm salinity model (qvalue <0.05)
* `LEA3.2_lfmm2_num_causal_sig_temp` number of causal loci on the temp trait, significant in the lfmm temp model (qvalue <0.05)
* `LEA3.2_lfmm2_num_neut_sig_temp` number of neutral loci false positives (only neutral loci not affected by selection), significant in the lfmm temp model (qvalue <0.05)
* `LEA3.2_lfmm2_num_causal_sig_sal` number of causal loci on the salinity trait, significant in the lfmm salinity model (qvalue <0.05)
* `LEA3.2_lfmm2_num_neut_sig_sal` number of neutral loci false positives (only neutral loci not affected by selection), significant in the lfmm salinity model (qvalue <0.05)
*  `LEA3.2_lfmm2_FPR_neutSNPs_temp` false positive rate of lfmm temperature model
*  `LEA3.2_lfmm2_FPR_neutSNPs_sal` false positive rate of lfmm salinity model
*  `RDA1_temp_cor` output of `summary(rdaout)$biplot[2,1]`, which is the correlation between RDA1 and the temperature environmental variable
*  `RDA1_sal_cor`  output of `summary(rdaout)$biplot[1,1]`, which is the correlation between RDA1 and the salinity environmental variable
*  `RDA2_temp_cor` output of `summary(rdaout)$biplot[2,2]`, which is the correlation between RDA2 and the temperature environmental variable
*  `RDA2_sal_cor` output of `summary(rdaout)$biplot[1,2]`, which is the correlation between RDA2 and the salinity environmental variable
*  `RDA_mlog10P_sig_noutliers` number of outliers in the RDA analysis 
*  `RDA_FPR_neutSNPs` false positive rate of the RDA analysis, based only on neutral SNPs
*  `RDA_RDA1Loading_cor_tempEffect` correlation between loading of a mutation on the first RDA axis (usually temperature) and effect size of allele on temperature
*  `RDA_RDA1Loading_cor_salEffect` <- correlation between loading of a mutation on the first RDA axis (usually temperature) and effect size of allele on salinity
*  `RDA_RDA2Loading_cor_tempEffect` <- correlation between loading of a mutation on the second RDA axis (usually salinity) and effect size of allele on temperature
*  `RDA_RDA2Loading_cor_salEffect` <- correlation between loading of a mutation on the second RDA axis (usually salinity) and effect size of allele on salinity
  
  something to think about - this RDA approach might work only because ancestral alleles are neutral.

Decided not to do
  - [ ] change "neutral-linked" to "non-causal"
  - [ ] add INFO to mutations output (only adds pop origin and gen origin)
  - [ ] need to add mutation-specific and genome-wide FST  calculation to output and outliers


- [ ] **Parameters**
  - [ ]  Previously I got cool results with the polygenic mutation rate with Sigma_QTN=0.1. Now I have it set to sigma_QTN=0.002. The oligogenic case is set to Sigma_QTN=0.4. So I think we should revise the parameter set 
  - [ ]  Do higher gene flow - less LA?


- [ ] **Housekeeping**
  - [ ] Download YML files from `/work/lotterhos/MVP-NonClinalAF/src` to  https://github.com/northeastern-rc/lotterhos_group

