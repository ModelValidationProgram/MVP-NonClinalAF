# MVP Repo

Files
* 0a-writeSimScript.R
* LFMM_tutorial.R
* README.md 

Folder
* `img/` - Image files
* `notebook` - Notebook entries
* `output` - Active outputs
* `output_archived` - Archived outputs 
* `slurm_log` - Slurm output logs
* `src` - Scripts
* `sim_output_150_20210830` - outputs from an intial set of SliM runs. Some of the runs may have bugs and some may not, as code was being updated while these were running. This had simulations with m_breaks set to 0.0001, which caused some of the estuary sims to take a really long time to coalesce.
* 

## Environments
MVP_env for running SLiM
MVP_env_R4.0.3 for running R

## Preparing a run
In `d-run_nonAF_sims_0Slim.sh` specify:
* `jobname`, `output`, and `error` files at the top
* `outpath` 
* `params` file 

Use a naming system with the date that corresponds to the notebook post date.

For low memory jobs:
```
#SBATCH --mem=5G
#SBATCH --nodes=1
#SBATCH --array=2-1000%72
```
Each core can have up to 5G memory, and an array won't run more jobs than cores available. So up to 72 jobs in the array can be requested at a time.

From /work/lotterhos/MVP_NonClinalAF/ run the SliM Script: `sbatch d-run_nonAF_sims_0Slim-fastruns-20220117.sh`

Check the run: `squeue -u lotterhos`

Check what's running on the lotterhos nodes: `squeue -p lotterhos`

Make sure the outputs are there and then edit the R script:

TO DO

Run the R script:

