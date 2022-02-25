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
In `d-run_nonAF_sims_0Slim.sh` copy it to a file name with the date, e.g. `d-run_nonAF_sims_0Slim-fastruns-20220117.sh` and specify:
* Use the same date in the `jobname`, `output`, and `error` files at the top
* Use the same date in the `outpath` 
* Use the same date in the `params` file 


Use a naming system with the date that corresponds to the notebook post date.

For low memory jobs:
```
#SBATCH --mem=5G
#SBATCH --nodes=1
#SBATCH --array=2-1000%70
```
Each core can have up to 5G memory, and an array won't run more jobs than cores available. So up to 72 jobs in the array can be requested at a time. I like to leave a couple open for `srun`.

From `/work/lotterhos/MVP_NonClinalAF/` run the SliM Script: `sbatch d-run_nonAF_sims_0Slim-fastruns-20220117.sh`

Check the run: `squeue -u lotterhos`

Make sure to write the jobID in the notebook post

Check what's running on the lotterhos nodes: `squeue -p lotterhos`

Make sure the outputs are there and then copy the R script `e-run_nonAF_sims_1R.sh` to a file name with the date, e.g. `e-run_nonAF_sims_1R-fastruns-20220117.sh` 
* Use the same date in the `jobname`, `output`, and `error` files at the top
* Use the same date in the `outpath` 
* Use the same date in the `params` file 

The date is like a runID, so using the same date across all files helps to coordinate different runs.

Edit the R script `c-AnalyzeSimOutput.R` to update with the date of the runs on this line:
`allsims <- load("src/0b-final_params-20220117.RData")`

Run the R script:

`sbatch e-run_nonAF_sims_1R-fastruns-20220117.sh`

Check the runs: `squeue -u lotterhos`

Make sure to write the jobID in the notebook post

To cancel the runs: `scancel <jobID>`

After the runs are finished, check how efficient they were with `seff <jobID>` and write that info in the notebook post

https://rc-docs.northeastern.edu/en/latest/using-discovery/usingslurm.html

## Git tips on Discovery

```
vi .gitignore #edit the gitignore file
git branch # make sure you are on the right branch!
git status # gives an overview of what's going on
git add -A # add/stage all changes to be commited
git commit # commit the changes
git pull # pull any updates from the remote repo (do not run this before you commit!)
git push # push the committed changes to the repo
```
