# Description of final shell script file - 2021-Jul-01



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
