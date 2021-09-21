

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

# 2021 09 20

## update parameter list

- [x] **Parameters**
  - [x]  Previously I got cool results with the polygenic mutation rate with Sigma_QTN=0.1. Now I have it set to sigma_QTN=0.002. The oligogenic case is set to Sigma_QTN=0.4. So I think we should revise the parameter set 
    - [x]  Sigma_QTN=0.4, proportion of mutations = 0.001 (same as before)
    - [x]  sigma_QTN=0.002, proportion of mutations = 0.25 (same as before)
    - [x]  Sigma_QTN=0.1, proportion of mutations = 0.1 (set to add)
  - [ ]  Do higher gene flow - less LA - over a larger number of sims?

I created a new parameters file, we now have 225 levels in the simulations, each replicated 10 times.

I separated the runs into "long runs" and "fast runs".

`longruns <- which(final$demog_name=="Est-Clines_N-variable_m-variable" | final$demog_name=="Est-Clines_N-equal_m_breaks")`

* `src/0b-final_params-fastruns-20210920.txt` Give these standard memory and time limits
* `src/0b-final_params-longruns-20210920.txt` Give these more memory and time limits

In addition, I saved the params file as an R object (since I had all the factors ordered appropriately when I made it)

`src/0b-final_params-20210920.RData`

# start runs

I discovered that the max array size is 1000! That's no fun. 

The fast runs ahve 1950 sims. The long runs have 300 sims.

Job IDS: 

21152261_2, 21152272 (array index 3 to 1000, first 1000 parameters),  21174128 (array index 2 to 1950, 1000-1950 parameter levels), 21174205 (long sims levels)


### sbatch: error: Batch job submission failed: Invalid job array specification

This error arose anytime I used an array index greater than 1000, for example:
`#SBATCH —array=1001-1951%68`  because the array index goes over 1000 (Even though the number of jobs in the array is less than 1000)

So, I split up the sims into 3 datasets, each submitted as a separate array. 

I also noticed that the lotterhos partition could probably handle more sims then was allowed in the array limit (68 jobs submitted at a time). I think there is 5GB of memory on each core, with 72 cores total across both nodes. So I think I could submit an array with 140 jobs at a time, each with 2.5 GB memory, but I'm not sure.

For 68 jobs of 2G allocated each I get:
```
$scontrol show node d3037

NodeName=d3037 Arch=x86_64 CoresPerSocket=18
   CPUAlloc=36 CPUTot=36 CPULoad=22.42
   RealMemory=190000 AllocMem=92160 FreeMem=159527 Sockets=2 Boards=1
   State=ALLOCATED ThreadsPerCore=1
   Partitions=lotterhos
   CfgTRES=cpu=36,mem=190000M,billing=36
   AllocTRES=cpu=36,mem=90G
```
CfgTRES=cpu=72,mem=1519000M,billing=6005
AllocTRES=cpu=72,mem=441840M

It shows that for this node, it has 36 cores and 190 GB (0.19TB) RAM (the CfgTRES part), and currently jobs are using 36 cores, and 90GB RAM (the AllocTRES part).

190 GB/ 36 cores = 5.27 GB per CPU.

So I could submit more jobs at a time with the array.

# To Do 

- [ ] when first round of short jobs finish, run second round of short jobs
- [ ] create array submission for long jobs

- [ ] **Housekeeping**
  - [ ] Download YML files from `/work/lotterhos/MVP-NonClinalAF/src` to  https://github.com/northeastern-rc/lotterhos_group

- [ ] **Methods**
  - [ ] update methods of manuscript

