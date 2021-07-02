# Description of final shell script file - 2021-Jul-01

## Overview 

Shell script has been updated to leverage slurm `--arrays` and read parameters from a [parameter file](https://github.com/ModelValidationProgram/MVP-NonClinalAF/blob/alan/src/0b-final_params.txt). It does this by temporary storing variables from the parameter file in unix variables, which are used by SLiM and other programs in the script. Below I provide a basic walk through of the different parts of the script.


**Final Shell Script** : [src/run_nonAF_sims_alan_final.sh](https://github.com/ModelValidationProgram/MVP-NonClinalAF/blob/alan/src/run_nonAF_sims_alan_final.sh)

### Test outputs

Test outputs can be found on the discovery cluster at:

```
/work/lotterhos/MVP-NonClinalAF/finalTest
```

For the test I ran the first 10 samples simultaneuosly on the short partition. To make the script run faster I limited the population size in pyslim to `10`. I checked that the vcf looked "OK" by reading them into R with `vcfR`.


### Slurm Code

Slurm code is standard, but it does use the `--array` argument. This can be used to run through a subset or all of the parameters combinations (rows) in your parameter list.

```
#SBATCH --job-name=SlimTest
#SBATCH --mem=50G
#SBATCH --mail-user=downey-wall.a@northeastern.edu
#SBATCH --mail-type=FAIL
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --array=2-5%4
#SBATCH --output=/work/lotterhos/MVP-NonClinalAF/slurm_log/SlimTest_%j.out
#SBATCH --error=/work/lotterhos/MVP-NonClinalAF/slurm_log/SlimTest_%j.err
```

**Couple of notes about the array**:

1) The array should ALWAYS start on 2. This is because your want to skip the first row of the parameter text files that contains the column names.
2) The largest value in the array range should not exceed the number of rows in your parameter list.
3) The number of simulations running simultaneously,`%4`, should not exceed the total number of rows you are looping through with the array.

**KEL question** So in the above script, you are going to run rows 2-5 in the params file (4 sims), each simultaneously?

### Custom user variables

At the top of the script (after the slurm arguments) I have the user inputs:

```
#### User modified values ####

# Local working path (this should navigate to the MVP repo)
mypath="/work/lotterhos/MVP-NonClinalAF"
cd ${mypath}
# Folder within MVP where you want are your output files
outpath="sim_outputs_testAlanV2/"
mkdir -p ${outpath} # make outpath directory if it doesn't exist

# Parameter file
params="src/0b-final_params.txt"

#### User variables ####
# N for pyslim
POP=10
# Minimum allele freq. for vcftools
MAF=0.01
```

Here you can modify the (i) path to the MVP-NonClinalAF repo as needed (`mypath`), (ii) the folder within the MVP-NonClinalAF repo where you want your outputs (`outpath`), and couple of pyslim and vcftools related parameters which were originally hardcoded in the script (`POP` and `MAF`).

Note on the `outpath` - The `outpath` will automatically be created if it does not already exist.

### Parameter variables

In this version of the shell script the SLiM and pyslim parameters are read directly from the parameter txt file and stored temporary as unix variables. These variables are used later in the script. The code chunk looks like this:

```
level=`awk NR==${SLURM_ARRAY_TASK_ID} ${params} | awk '{print $1}'`
reps=`awk NR==${SLURM_ARRAY_TASK_ID} ${params} | awk '{print $2}'`
arch=`awk NR==${SLURM_ARRAY_TASK_ID} ${params} | awk '{print $3}'`
demog_name=`awk NR==${SLURM_ARRAY_TASK_ID} ${params} | awk '{print $4}'`
demog_level_sub=`awk NR==${SLURM_ARRAY_TASK_ID} ${params} | awk '{print $5}'`
demog_level=`awk NR==${SLURM_ARRAY_TASK_ID} ${params} | awk '{print $6}'`
MIG_x=`awk NR==${SLURM_ARRAY_TASK_ID} ${params} | awk '{print $7}'`
MIG_y=`awk NR==${SLURM_ARRAY_TASK_ID} ${params} | awk '{print $8}'`
xcline=`awk NR==${SLURM_ARRAY_TASK_ID} ${params} | awk '{print $9}'`
ycline=`awk NR==${SLURM_ARRAY_TASK_ID} ${params} | awk '{print $10}'`
demog=`awk NR==${SLURM_ARRAY_TASK_ID} ${params} | awk '{print $11}'`
METAPOP_SIDE_x=`awk NR==${SLURM_ARRAY_TASK_ID} ${params} | awk '{print $12}'`
METAPOP_SIDE_y=`awk NR==${SLURM_ARRAY_TASK_ID} ${params} | awk '{print $13}'`
Nequal=`awk NR==${SLURM_ARRAY_TASK_ID} ${params} | awk '{print $14}'`
isVariableM=`awk NR==${SLURM_ARRAY_TASK_ID} ${params} | awk '{print $15}'`
MIG_breaks=`awk NR==${SLURM_ARRAY_TASK_ID} ${params} | awk '{print $16}'`
arch_level_sub=`awk NR==${SLURM_ARRAY_TASK_ID} ${params} | awk '{print $17}'`
arch_level=`awk NR==${SLURM_ARRAY_TASK_ID} ${params} | awk '{print $18}'`
MU_base=`awk NR==${SLURM_ARRAY_TASK_ID} ${params} | awk '{print $19}'`
MU_QTL_proportion=`awk NR==${SLURM_ARRAY_TASK_ID} ${params} | awk '{print $20}'`
SIGMA_QTN_1=`awk NR==${SLURM_ARRAY_TASK_ID} ${params} | awk '{print $21}'`
SIGMA_QTN_2=`awk NR==${SLURM_ARRAY_TASK_ID} ${params} | awk '{print $22}'`
SIGMA_K_1=`awk NR==${SLURM_ARRAY_TASK_ID} ${params} | awk '{print $23}'`
SIGMA_K_2=`awk NR==${SLURM_ARRAY_TASK_ID} ${params} | awk '{print $24}'`
N_traits=`awk NR==${SLURM_ARRAY_TASK_ID} ${params} | awk '{print $25}'`
ispleiotropy=`awk NR==${SLURM_ARRAY_TASK_ID} ${params} | awk '{print $26}'`
seed=`awk NR==${SLURM_ARRAY_TASK_ID} ${params} | awk '{print $27}'`
```

**NOTE** - If you change the parameters in your parameters (and SLiM model) you can easily update this block of code by walking through the steps from [this notebook entry](https://github.com/ModelValidationProgram/MVP-NonClinalAF/blob/alan/notebook/20210630_creatingSLiMBatchScriptVariables.md).

### SLiM code

This didn't really change, but I did switch the hardcoded values with the temporary unix variables that represent values from the parameter list.

```
slim -d MY_SEED1=${seed} -d MY_RESULTS_PATH1=\"${outpath}\" -d MIG_x1=${MIG_x} -d MIG_y1=${MIG_y} -d demog1=${demog} -d xcline1=${xcline} -d ycline1=${ycline} -d METAPOP_SIDE_x1=${METAPOP_SIDE_x} -d METAPOP_SIDE_y1=${METAPOP_SIDE_y} -d Nequal1=${Nequal} -d isVariableM1=${isVariableM} -d MIG_breaks1=${MIG_breaks} -d MU_base1=${MU_base} -d MU_QTL_proportion1=${MU_QTL_proportion} -d SIGMA_QTN_1a=${SIGMA_QTN_1} -d SIGMA_QTN_2a=${SIGMA_QTN_2} -d SIGMA_K_1a=${SIGMA_K_1} -d SIGMA_K_2a=${SIGMA_K_2} -d N_traits1=${N_traits} -d ispleiotropy1=${ispleiotropy} src/a-PleiotropyDemog.slim > sim_outputs_testAlan/${seed}_outfile.txt 2> sim_outputs_testAlan/${seed}_outfile.error.txt    
```

### pyslim

Similar to SLiM code, I didn't change much here but replaced hardcode with temporary variables. I also moved the population size variable, `N`, up to the top of the script. It is now controlled by a user variable called `POP`.

```
python3 src/b-process_trees.py -s ${seed} -r 1e-06 -mu ${MU_base} -N ${POP} -o ${outpath} > ${outpath}${seed}_pytree.out.txt 2> ${outpath}${seed}_pytree.error.txt
```

### VCF files and filtering

Only minor updates here. I moved the minimum allele frequency value up to the top of the script. It is not controlled by a user variable called `MAF`.

```
vcftools --vcf ${outpath}${seed}"_PlusNeuts".vcf --maf ${MAF} --out ${outpath}${seed}"_plusneut_MAF01" --recode --recode-INFO-all
```

### Timestamps

I slightly revised the timestamp system. These shouldn't need to be changed, but leverage the built in `date` function in unix to create a timestamp and print how long the code had been running. This code generally follows the format:

**Start timestamp**
```
echo "Running SLiM..."
timestamp_initial=`date`
echo ${timestamp_initial}
```

**End timestampe and print**
```
## Final timestamp
echo "Done with SLiM sims"
time_stampFinishSlim=`date`
echo ${time_stampFinishSlim}
HOURS=$(echo $(date -d "$time_stampFinishSlim" +%H) - $(date -d "$timestamp_initial" +%H) | bc)
MINS=$(echo $(date -d "$time_stampFinishSlim" +%M) - $(date -d "$timestamp_initial" +%M) | bc)
SECS=$(echo $(date -d "$time_stampFinishSlim" +%S) - $(date -d "$timestamp_initial" +%S) | bc)
echo "Slim run took: ${HOURS} hours, ${MINS} minutes, and ${SECS} seconds"
```

It updates the user after SLiM, pyslim, and the script completes.



