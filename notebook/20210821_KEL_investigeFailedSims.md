Alan ran the first 50 simulations in the parameter file, and a number of them didn't finish within 2 weeks. 
I hypothesize that this issue arose due to failure of the pyslim script to coalesce, which has happened in the past when the migration matrix was not
properly specified.

In the:

20210818_150_V2_SummaryReport.html put together by Alan, all of the demographies that failed were Estuary.

On Discovery I looked at the outputs in:

`/work/lotterhos/MVP-NonClinalAF/sim_output_150_V2`

`ls -lah`
revealed an interesting pattern - the dates the files are written seem to be in succession. Why is that? If the partition wasn’t busy, shouldn’t they have all started at the same time? (maybe this is why some of them timed out - the 2 weeks is from the time the array script was submitted)

Thoughts:
- 1. I need to check if the VCF files that I'm producing are too large to work with, even after filtering. Maybe reduce the population size or mutation rate in pyslim. File size increases from 30-40KB to 322MB
  - try opening VCF files in VCFR - there are some edits I wanted to make anyway
- 2. check that pyslim can work on the failed sims to produce a VCF file.
  - 1231143, 1231141, 1231096, 1231098 as starting points

## Summary of failed sims
The failed simulations include:
* 5/5 Estuary-Clines N-equal_m_breaks
* 4/5 Estuary-Clines N-variable_m-variable

The fact that one of the N-variable_m-variable finished makes me think it is taking a long time for the sims to coalesce, and it's not an issue with the sims.


# 2021-08-23

## Setup
```
ssh lotterhos@discovery.neu.edu
cd /work/lotterhos/MVP-NonClinalAF/

ssh -T git@github.com # not sure if I needed this
git branch # just checking

srun -p lotterhos -N 1 --pty /bin/bash

source ~/miniconda3/bin/activate MVP_env # activate the environment
python3
```
In python:
```
import pyslim 
import msprime
import argparse
import numpy as np
import random
import time
import re
import ast
import sys
import os
seed = '1231143'
r = 1e-6
N = 10
mu = 1e-7
gen = 3000
seednum = round(int(seed)**(1/2)) # slim seeds are too large
output = ""
os.chdir("/work/lotterhos/MVP-NonClinalAF/sim_output_150_V2")
T2 = pyslim.load(output + seed + ".trees")
```
This gives me an error: `Tree sequence contains no SLiM provenance entries(or your pyslim is out of date).`

Back in command line:
```
conda update pyslim
python3
```
I ran the above python lines, but still the error: `Tree sequence contains no SLiM provenance entries(or your pyslim is out of date).`


Alan: "as a backup I made a new environment text file from my current MVP_env conda environment. You should be able to create a new MVP_env that mirrors my current environment with the following command:"
```
cd /work/lotterhos/MVP-NonClinalAF/src/env
conda create --name MVP_env --file MVP_env_V2.tx
```
I overwrote the old MVP_env with this new environment

Now the python lines are working!

## Testing recapitation on sims that didn't finish

I tested recaptiation with pyslim for `N=1` and `r=1e-11` for a sim that didn't finish for Alan (1231143) and it worked! Therefore, there's not a problem with the sims, but with the time limits on the cluster.

## Testing if vcf files work on my R code in OOD

Alan ran the first 50 sims, which did not include any of the polygenic sims or sims with pleiotropy.

I launched R Studio, R 4.0.3 with 32 GB of memory (2x my laptop) in OOD

Testing on `seed = 1231107`, oliogenic_1-trait__SS-Mtn_N-equal_m-constant

**The VCF files produced after filtering for MAF are 350MB large, which is kind of big. When I load it into R, it needs 10Gb of memory**

Started reading in file to R studio 10:38am. By 10:45 it starts to read in the loci. By 10:51 it's halfway through loading the loci. Some of the slowness comes from the large number of 10000 individuals, but it has to happen this way so I can sample according to fitness.

  variant count: 133084
  column count: 10009

OOD annoyance: After not being active for a bit, my Rstudio session needed to reload. However it had that huge VCF file that took several minutes to load, so it was pretty unhappy.cp


## Interesting observation about git on Discovery

When I'm in the login node, git seems to work. 
```
ssh -T git@github.com
Hi ModelValidationProgram/MVP-NonClinalAF! You've successfully authenticated, but GitHub does not provide shell access.
```

**When I'm in an srun environment, the above line of code does not seem to work.
This means that any changes to git files have to happen in the login node.**

OR: I need to redo the connection for the lotterhos node.

IN ADDITION: When I try to push, it doesn't work. It says: "Support for password authentication was removed on August 13, 2021. Please use a personal access token instead.
remote: Please see https://github.blog/2020-12-15-token-authentication-requirements-for-git-operations/ for more information."

I followed the steps on this page: https://docs.github.com/en/github/authenticating-to-github/keeping-your-account-and-data-secure/creating-a-personal-access-token

I stored the token on my computer. Now when I push to git, I use this as a password.

## Interesting observation about png in R on Discovery

For some reason I get an error for png, but not pdfs. ignoring for now, just outputting pdfs.

## RDA note

I increased the number of loci in the RDA that predicts the traits to 50,000. That may be too many if I decrease the N from 1000 to 100 for pyslim - we'll see!

## R version and Conda notes

In debugging my R code on OOD, I was using R version 4.0.3
However, the conda environment only has R version 3.6

Here is how I updated the environment:
```
conda create -n MVP_env_R4.0.3
conda config --add channels conda-forge
conda config --set channel_priority strict
conda search r-base
conda activate MVP_env_R4.0.3
conda install -c conda-forge r-base=4.0.3
conda list
conda env export -f MVP_env_R4.0.3.yml
```

environment location: /home/lotterhos/miniconda3/envs/MVP_env_R4.0.3


After activating this environment, I had to install all the libraries that I needed. I typed `R` at the command line, installed the packages, and check the
library location with `.libPaths()`, which returned:
`/home/lotterhos/miniconda3/envs/MVP_env_R4.0.3/lib/R/library`

Just to be safe, after I intalled the R packages I quit R and I exported the environment again with:
`conda env export -f MVP_env_R4.0.3.yml`

YML location: `/work/lotterhos/MVP-NonClinalAF/src/env/MVP_env_R4.0.3.yml`

## TO DO

-[x] mountain range sims suggest temperature at demes 81-90 is higher than 91-100. FIX SLIM. (Fixed - see below)

-[x] make an R script based on the markdown
  - [x] OOD .libPaths() "/home/lotterhos/R/x86_64-pc-linux-gnu-library/4.0" "/shared/centos7/r-project/4.0.3/lib64/R/library"
  - [x] CONDA "/home/lotterhos/miniconda3/envs/MVP_env_R4.0.3/lib/R/library"
  
-[x] incorporate R code into bash script


## TO DO

-[x] pyslim recaptiation is a slow step. Alan confirmed that he did run it with N=1000. It also produces huge VCF files.

R Studio on chunck  - I left off on chunck 27 with heatmaps


-[x] New Submission script: `run_nonAF_sims_0Slim_20210826.sh`
  - Renamed to date
  - Renamed outputs to date
  - Used N=100 for pyslim
  - Run 2-151, we have 70 cores on Lotterhos

## First submission
- [x] Submit script for 150 sims
  - [x] `sbatch run_nonAF_sims_0Slim_20210826.sh`
    - trying 2 nodes, 70 per node 
    - should I specify --cpus-per-task=1?
    - For some reason it only started running the first 6 jobs, which was annoying
  - [x] check that it started `squeue -u lotterhos`
  - [x] check all jobs on lotterhos partiition with `squeue -p lotterhos`
  - Job ID: Submitted batch job 20695234
  - [x] check job efficiency `seff 20695234`
    - `badly formatted array jobid 20695234 with task_id = -2` Probably not enough memory
    -  that's not good!
  -  check specifics `scontrol show jobid -dd 20714872`
      -   `NumNodes=2-2 NumCPUs=2 NumTasks=2 CPUs/Task=1 ReqB:S:C:T=0:0:*:*
   TRES=cpu=2,mem=60G,node=2,billing=2`

- [x] Email Research Computing
- [x] Meet with RC
  - [x] Want to upload R to version 4.0.3 in my conda environment /work/lotterhos/MVP-NonClinalAF/src/env/MVP_env.yml or make a new environment
  - [x] Want to make submission of simulations more efficient: /work/lotterhos/MVP-NonClinalAF/src/run_nonAF_sims_0Slim.sh 
  - [ ] How to set the maximum time limits (if needed)
  - [x] submit R script with dependency that previous script finishes first


### 8/27 second submission
- the first 48 sims ran, but the same 6 simulations that took a long time before are also taking a long time
- I canceled the job `scancel 20695234` so I could adjust the settings
- I edited `run_nonAF_sims_0Slim.sh` and removed the memory allocation
- I submitted `sbatch run_nonAF_sims_0Slim_20180827.sh`
  - job ID `20714872`  
  - still writing to folder `sim_output_150_20210826`
- going to start sims from line 50 (sim 49)
- already this is looking better, 34 runs started at once
- `seff 20714872` still get `badly formatted array message`
- `scontrol show jobid -dd 20714872`
    - `NumNodes=2 NumCPUs=2 NumTasks=2 CPUs/Task=1 ReqB:S:C:T=0:0:*:*
   TRES=cpu=2,mem=4000M,node=2,billing=2`

### 8/27 initial results from sims with pyslim recap N=100

- ~20K as opposed to 300K+ variants (maybe N=500 would be a nice compromise)
  - 17Mb largest object  
- stepping stone has issues with opt1. I need to debug these! the -1 is missing somewhere.

### SLIM ISSUES FIXED

- There was an issue with setting the optimums in the slim simulation. I fixed it and pushed to github. This means that the simulations in the output folder will all have slightly incorrect optimums: `sim_output_150_20210826`, but will all be good for de-bugging. FIXED AND SYNCED WITH GIT.

- Also a bug in the setting the genomic element - I had set it to the first 11 elements (LGs), but only mean to set it to the first 10. FIXED AND SYNCED WITH GIT.

## R script edits
- [x] For the three methods, OUTPUT for outliers: VA_temp, VA_sal, FDR_allSNPs, FDR_neutSNPs, AUCPR_allSNPs, AUCPR_neutSNPs, and manhattan plots
  - [x] correlation
  - [x] RDA
  - [x] LFMM
- [x] ADD COR(AF, PC1) to understand effect of correction for sturcture
- [x] ADD plot to understand if RDA loading corresponds to / locus effect size

## Test R Bash Script
- [x] Make an R Bash script `run_nonAF_sims_1R.sh`
- [ ] Get environment with R v 4.0.3. Edit and submit for one simulation, see if it works!


## New Bug
  - [X] For 2 trait architucture, subsampling does not result in AF clines from Slim Correlating with AF clines from subsample for Env2, but it does for temperature. FIXED Aug 30
    - [x] BUG: when merging the indPhen_df, individual orders were mixed up. FIXED. This was a major bug affecting many outputs.



## Notes from RC meeting Aug 30



`#SBATCH --array=50-151%70`

* If each task in array requires 1 CPU, %72 to maximize lotterho partition.
* If it requires 2 CPUS, set --CPUs-per-task=2 (look it up) and set --array=50-151%36 
* IF programs can use multiple threads, do benchmarking to determine CPU per task, 1-2 jobs at a time, and do `seff` on those tasks that completed, test 1, 2, 4, 8, 16 CPUs
* Need to know number of CPUS per task, For maximum efficiency 



`#SBATCH --nodes=2`

* Set this to 1 unless we know the program can communicate between the nodes.


```
seff 20714872
Job ID: 20714872
Array Job ID: 20714872_151
Cluster: discovery
User/Group: lotterhos/users
State: COMPLETED (exit code 0)
Nodes: 2
Cores per node: 1
CPU Utilized: 00:13:00
CPU Efficiency: 48.51% of 00:26:48 core-walltime
Job Wall-clock time: 00:13:24
Memory Utilized: 768.23 MB
Memory Efficiency: 19.21% of 3.91 GB
```

`#SBATCH --mem=170GB`

This report is for each array task.
Each requires 1 GB memory.
So for the total array submission, mulitply mem/task * 

If I set `#SBATCH --array=50-151%70`, and each task needs 1GB memory, I would specificy total mem = 70*1GB = 70GB, but add buffer and set it to mem=170GB in case one task requires more memory.

If it doesn't specify anthing, it will use the default amount, for one task it might work but it might not work every time. The default memory is a limit, so it's always recommended to set it higher so you don't run out.ls

## Aug 30 new bash submission
- I synced Github with cluster. Hopefully I didn't F anything up.
- Still working on making sure R outputs are OK
- I canceled the last job `20714872`, since then I have fixed some of the slim scripts and some of the R outputs. I also learned how to improve the submission script from my meeting with RC
- created new submission for 20210830
- `sbatch run_nonAF_sims_0Slim_20210830.sh`
- Job ID 20883113
- `squeue -u lotterhos`
- `seff 20883113`
- `scontrol show jobid -dd 20883113`

Note: if I set mem=170G, it only runs one job. If I take out the memory option, it uses what it can of the cluster. If I set mem=2G, it uses what it can of the lotterhos partition. I think each CPU has memory up to 3.91 GB (based on seff output).
    
## Done
- [x] Rerun first 150 sims with revised SliM Code (fixed issues with optimum and genomic element)
- [x] R code need to remove rare neutral alleles <0.01 after sampling
- [x] **Add to outputs**
  - [x]  numCausalLowMAFsample
  - [x]   Bonf_alpha
  - [x]  subsamp_corr_phen_temp, subsamp_corr_phen_sal
  - [x]  K
  - [x]      cor_PC1_temp <- cor(subset_indPhen_df$PC1, subset_indPhen_df$temp_opt)
    cor_PC1_sal <- cor(subset_indPhen_df$PC1, subset_indPhen_df$sal_opt)
    cor_PC2_temp <- cor(subset_indPhen_df$PC2, subset_indPhen_df$temp_opt)
    cor_PC2_sal <- cor(subset_indPhen_df$PC2, subset_indPhen_df$sal_opt)
    cor_LFMMU1_temp <- cor(subset_indPhen_df$LFMM_U1_temp, subset_indPhen_df$temp_opt)
    cor_LFMMU1_sal <- cor(subset_indPhen_df$LFMM_U1_sal, subset_indPhen_df$sal_opt)
    cor_LFMMU2_temp <- cor(subset_indPhen_df$LFMM_U2_temp, subset_indPhen_df$temp_opt)
    cor_LFMMU2_sal <- cor(subset_indPhen_df$LFMM_U2_sal, subset_indPhen_df$sal_opt)
    cor_PC1_LFMMU1_temp <- cor(subset_indPhen_df$LFMM_U1_temp, subset_indPhen_df$PC1) 
    cor_PC1_LFMMU1_sal <- cor(subset_indPhen_df$LFMM_U1_sal, subset_indPhen_df$PC1) 
    cor_PC2_LFMMU1_temp <- cor(subset_indPhen_df$LFMM_U1_temp, subset_indPhen_df$PC2) 
    cor_PC2_LFMMU1_sal <- cor(subset_indPhen_df$LFMM_U1_sal, subset_indPhen_df$PC2) 

# Testing R script and R v4.0.3 Conda environment
For some reason, the path that works for `source ~/miniconda/bin/activate MVP_env`
doesn't exist anymore, the folder is now `miniconda3`

By using `conda info --envs` I found the locations of the environments are at:
`source ~/miniconda3/bin/activate MVP_env_R4.0.3`

`sbatch run_nonAF_sims_1R_20210830.sh`

20891462

`squeue -u lotterhos`

debugging:
```
vi ../slurm_log/R-Run20210830_20888XXX.err
vi /work/lotterhos/MVP-NonClinalAF/sim_output_150_20210830/1231094_R.error
vi /work/lotterhos/MVP-NonClinalAF/sim_output_150_20210830/1231094_R.out
```

Performance:
`seff 20891462_2`
```
Cores: 1
CPU Utilized: 00:19:20
CPU Efficiency: 99.83% of 00:19:22 core-walltime
Job Wall-clock time: 00:19:22
Memory Utilized: 8.71 GB
Memory Efficiency: 435.70% of 2.00 GB
```

# Running full R array

`sbatch run_nonAF_sims_1R_20210830.sh`

20891805


## figure out why these sims take so long: `Est-Clines_N-equal_m_breaks` `Est-Clines_N-variable_m-variable`

One of them finished! But it took longer than a week.

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

- [ ] 
  - [ ] recaptitate these with N=1 and r=1e-11, it won't take long to finish. Go in R to see what's going on.   
  - OR: Just report that they do not coalesce, and remove from the results
  - OR: don't let m get so small in the estuary demog and see if they recapitate
  - OR: it could be the default memory wasn't enough (Based on `seff` output I don't think this is the case

## To Do
- [x] R output for talk - compare 
  - [x] oligo SS 1 trait N-m equal to : 1231102
  - [x] oligo 2 traits plieotropy: 1231147 - 1231151
  - [x] polygenic SS 2 traits pleiotropy: 1231219 - 1231223
- [ ] **metadata for output dataframe**
- [ ] **Parameters**
  - [ ]  Previously I got cool results with the polygenic mutation rate with Sigma_QTN=0.1. Now I have it set to sigma_QTN=0.002. The oligogenic case is set to Sigma_QTN=0.4. So I think we should revise the parameter set 
- [ ] **R code** 
  - [ ] 1231150 gave an error in the histogram of the effect sizes - breaks do not span range of x
  - [ ] plot VA-prop vs. af cline
  - [ ] prop of causal alleles with significant clines?
  - [ ] change background for mutation AF clines density plot and explore different histogram types
  - [ ] need to add mutation-specific and genome-wide FST  calculation to output and outliers
- [ ] **Issues**
- [ ] Run R analysis for first 150 sims that worked
- [ ] Download YML files from `/work/lotterhos/MVP-NonClinalAF/src` to  https://github.com/northeastern-rc/lotterhos_group


