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
