# 20210623_InitialNotebookEntry_Alan

Today I worked on setting up the MVP repo on the Discovery:

### Updates

Changes to the repo:

1) Cloned MVP repo onto `/work/lotterhos`. For some reason I could not access the original version. I placed the old version temporary in a folder labeled `OLD`, in case there were things on there that are important
2) Incorporated some new folder into the MVP repo. Notably, I creates a `slurm_log` folder for all `.err` and `.out` slurm files when running sbatch. This is in the `.gitignore` folder, so it won't be on the remote version of the repo (although I can change this if we want this information accessible to everyone).


Changes to [`run_nonAF_sims_alan.sh`](https://github.com/ModelValidationProgram/MVP-NonClinalAF/blob/alan/src/run_nonAF_sims_alan.sh) script: 

1) Converted it into an sbatch script, with `#SBATCH` header arguments for designating specific resources. This will be expanded upon later to take advantages of arrays.
2) Created new CONDA env, `MVP_env`, which is sourced in the this script. This environment contains slim, python3.7, and vcftools (for the moment).
3) Cleaned up some of the time stamp code.
4) Begun replacing some of the hard coded paths with variables to improve flexibility. For example, "sim_output/" replaced with ${outpath}, where outpath="sim_output/"

Changes to [`b-process_trees.py`](https://github.com/ModelValidationProgram/MVP-NonClinalAF/blob/alan/src/b-process_trees.py) script:

1) I also added an argument to the script, `-o`, which allows the user to provide an output pathway to `.tree` script.
 
### Errors

Running this script I am getting an error when the python script is called:
```
ValueError: Tree sequence contains no SLiM provenance entries(or your pyslim is out of date).
```

This appears to be due to this line in the pyslim script:
```
T2 = pyslim.load(output + seed + ".trees")
```

I am not familiar enough with pyslim yet to know why this error is occuring.

### Note on CONDA env
To create a copy of the conda environment navigate to the `env` folder:
```
cd /work/lotterhos/MVP-NonClinalAF/src/env
```
and run 
```
conda env create -f MVP_env.yml
```
to create a copy of the environment.
