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



also request OOD for evaluating R code





