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
conda config --add channels conda-forge
conda config --set channel_priority strict
conda search r-base
conda install -c conda-forge r-base=4.0.3
```
It's taking a long time., then it didn't work when I typed `conda list` it still said R version was 3.6.0. Submit ticket to RC.


## TO DO

-[ ] mountain range sims suggest temperature at demes 81-90 is higher than 91-100. FIX SLIM. (I CAN'T FIND THE ISSU HERE, SEE IF PERSISTS IN OTHER CLINAL SIMS)

-[x] make an R script based on the markdown
  - [ ] .libPaths() "/home/lotterhos/R/x86_64-pc-linux-gnu-library/4.0" "/shared/centos7/r-project/4.0.3/lib64/R/library"
  
-[ ] incorporate R code into bash script


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
    - `badly formatted array jobid 20695234 with task_id = -2`
    -  that's not good!
  -  check specifics `scontrol show jobid -dd 20695234`
      -   `NumNodes=2-2 NumCPUs=2 NumTasks=2 CPUs/Task=1 ReqB:S:C:T=0:0:*:*
   TRES=cpu=2,mem=60G,node=2,billing=2`

- [x] Email Research Computing
- [ ] Meet with RC
  - Want to upload R to version 4.0.3 in my conda environment /work/lotterhos/MVP-NonClinalAF/src/env/MVP_env.yml
  - Want to make submission of simulations more efficient: /work/lotterhos/MVP-NonClinalAF/src/run_nonAF_sims_0Slim.sh 
  - How to set the maximum time limits
  - 


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
- stepping stone has issues with opt1. I need to debug these! the -1 is missing somewhere.

### SLIM ISSUE FIXED

- There was an issue with setting the optimums in the slim simulation. I fixed it and pushed to github. This means that the simulations in the output folder will all have slightly incorrect optimums: `sim_output_150_20210826`, but will all be good for de-bugging.


