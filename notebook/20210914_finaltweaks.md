

## Continuing to tweak `Est-Clines_N-equal_m_breaks` `Est-Clines_N-variable_m-variable`

These were sims that were producing a lot of SNPs and using a lot of memory. 
Low migration rates on the landscape resulted in the sims taking a long time to coalesce, 
and in addition a very large number of SNPs would evolve, which took a long time to filter and process.

I changed the migration rates for these scenarios in SliM, in the file: `d-run_nonAF_sims_0Slim_20210914.sh`, 
so there woudn't be quite so much isolation between populations.

```
(base) [lotterhos@login-01 ~]$ seff  21072516
Job ID: 21072516
Array Job ID: 21072516_141
Cluster: discovery
User/Group: lotterhos/users
State: COMPLETED (exit code 0)
Cores: 1
CPU Utilized: 01:34:18
CPU Efficiency: 97.87% of 01:36:21 core-walltime
Job Wall-clock time: 01:36:21
Memory Utilized: 1.46 GB
Memory Efficiency: 14.59% of 10.00 GB
```

Sims easier to handle! Checking R script...

## To Do 

- [ ] **Parameters**
  - [ ]  Previously I got cool results with the polygenic mutation rate with Sigma_QTN=0.1. Now I have it set to sigma_QTN=0.002. The oligogenic case is set to Sigma_QTN=0.4. So I think we should revise the parameter set 
    - [ ]  Sigma_QTN=0.4, proportion of mutations = 0.001 (same as before)
    - [ ]  sigma_QTN=0.002, proportion of mutations = 0.25 (same as before)
    - [ ]  Sigma_QTN=0.1, proportion of mutations = 0.1 (set to add)
  - [ ]  Do higher gene flow - less LA - over a larger number of sims?


- [ ] **Housekeeping**
  - [ ] Download YML files from `/work/lotterhos/MVP-NonClinalAF/src` to  https://github.com/northeastern-rc/lotterhos_group

- [ ] **Methods**
  - [ ] update methods of manuscript

