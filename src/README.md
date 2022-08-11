This is a README file for the `/src` directory.

The scripts are run in the order of their filenames.

* `0a-setUpSims.Rmd` A markdown file used to create the list of simulations to run. Outputs:
	* `0b-final_params-20220428.RData` and `0b-final_params-20220428.txt` are equivalent files, the former as an R object
	* `0b-final_params-20220428.txt_metadata.md` metadata for the above two files
	* related files only for dividing up runs on the cluster

* `a-PleiotropyDemog_20220428.slim` The SLiM code for all simulations.
* `b-process_trees.py` The pyslim code for recapitation and adding neutral mutations.
* `c-AnalyzeSimOutput.R` The R code that samples individuals from the metapopulation and conducts the pop gen analyses (clines, pca, lfmm, rda, etc)
* `d-run_nonAF_sims_0Slim.sh` The bash script that takes input from file `0a` and calls files `a` and `b`, then conducts VCF filtering.
	* related files only for dividing up runs on the cluster
* `e-run_nonAF_sims_1R.sh` The bash script that runs file `c` after `d` is finished

* `env` the conda environments for the runs

* `extraRcode.R` extra code not used for analysis
* `f-catoutputs.sh` cat together outputs from R runs
* `g-CheckRunsAreComplete.Rmd` extra code not used for analysis, was used for earlier debugging
* `g-FinalAnalysis.Rmd` the code used to create the plots for publication
* `parameterOrderForShellScript.txt` parameter order based on columns in `0b-final_params-20220428.txt`, beware hard coding
* `plot_mig_Graphs.R` some R code for plotting migration graphs for manuscript



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
