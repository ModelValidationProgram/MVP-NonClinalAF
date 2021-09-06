# Aug 2021

### To Do: figure out why these sims take so long: `Est-Clines_N-equal_m_breaks` `Est-Clines_N-variable_m-variable`

One of them finished! But it took a day. Others are still running.

Est-Clines_N-variable_m-variable

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

## To Do
- [ ] While I'm waiting for the slow sims to finish:
  - [x] **metadata for output dataframe** https://github.com/ModelValidationProgram/MVP-NonClinalAF/blob/alan/sim_outputs_meta/README.md
  - [ ] sync github
  - [ ] check outputs for other parameter levels and write a notebook post (e.g. m breaks scenarios)
  - [ ] Work on methods for manuscript
  - [ ] **Develop full analysis script**
    - [x] List of questions
    - [x] Brainstorm figures to answer the questions
    - [ ] cat together all the outputs
    - [ ] write R code to analyze and make figures
- [ ] **R code** 
  - [ ] 1231150 gave an error in the histogram of the effect sizes - breaks do not span range of x
  - [ ] change "neutral-linked" to "non-causal"
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

