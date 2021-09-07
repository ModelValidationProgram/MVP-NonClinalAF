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

- [ ] check how long 20883113_141 took
  - [ ] recaptitate these with N=1 and r=1e-11, it won't take long to finish. Go in R to see what's going on.   
  - OR: Just report that they do not coalesce, and remove from the results
  - OR: don't let m get so small in the estuary demog and see if they recapitate
  - OR: it could be the default memory wasn't enough (Based on `seff` output I don't think this is the case

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

#### Create new slim file a-PleiotropyDemog_202109.slim

The issue seems to be mostly with the Est-Clines_N-variable_m-variable demography, which needed a very large amount of memory and took a long time to run.

Migration rates were sampled from `sample(c(rep(0.001,6),0.001,0.01,0.1,0.25), 1)`, which I changed to `sample(c(0.001,0.01,0.1,0.25))`

I kept the sampling for the SS to `sample(c(rep(0.001,6),0.01,0.1,0.25)`, since this seemed to create cool patterns


## To Do

- [ ] While I'm waiting for the slow sims to finish:
  - [x] **metadata for output dataframe** https://github.com/ModelValidationProgram/MVP-NonClinalAF/blob/alan/sim_outputs_meta/README.md
  - [x] sync github
  - [x] check outputs for other parameter levels and write a notebook post (e.g. m breaks scenarios)
  - [ ] Work on methods for manuscript
  - [ ] **Develop full analysis script**
    - [x] List of questions
    - [x] Brainstorm figures to answer the questions
    - [ ] cat together all the outputs
    - [ ] write R code to analyze and make figures


- [ ] **R code** 
  - [ ] Errors
    - [ ] `Error in hist.default(muts_full$mutTempEffect, xlim = c(-1, 1), xlab = "Mutation temp effect",  : 
  some 'x' not counted; maybe 'breaks' do not span range of 'x'
Calls: hist -> hist.default`
    - [ ] `Warning message:
position_dodge requires non-overlapping x intervals 
Warning message:
position_dodge requires non-overlapping x intervals 
Error in mod_temp@U[, 2] : subscript out of bounds`

  - [ ] change "neutral-linked" to "non-causal"
  - [ ] export RDA pc loadings for individuals. I bet this is interesting for the m-variable case.
  - [ ] set salinity bubbles to be dark grey, decrease size of graphic output
  - [ ] RDA - add TPR for causal temp and causal sal
  - [ ] should I include the low MAF causal loci in the analysis? what about the highly polygenic case?
  - [ ] add INFO to mutations output
  - [ ] plot VA-prop vs. af cline
  - [ ] prop of causal alleles with significant clines?
  - [ ] change background for mutation AF clines density plot and explore different histogram types
  - [ ] need to add mutation-specific and genome-wide FST  calculation to output and outliers


- [ ] **Parameters**
  - [ ]  Previously I got cool results with the polygenic mutation rate with Sigma_QTN=0.1. Now I have it set to sigma_QTN=0.002. The oligogenic case is set to Sigma_QTN=0.4. So I think we should revise the parameter set 


- [ ] **Housekeeping**
  - [ ] Download YML files from `/work/lotterhos/MVP-NonClinalAF/src` to  https://github.com/northeastern-rc/lotterhos_group

