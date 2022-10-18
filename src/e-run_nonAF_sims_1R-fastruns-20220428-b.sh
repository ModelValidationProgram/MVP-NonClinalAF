#!/bin/bash

#SBATCH --job-name=R-Run-20221018
#SBATCH --mail-user=k.lotterhos@northeastern.edu
#SBATCH --mail-type=FAIL
#SBATCH --partition=lotterhos
#SBATCH --mem=15G
#SBATCH --nodes=1
#SBATCH --array=2-952%30
#SBATCH --output=/work/lotterhos/MVP-NonClinalAF/slurm_log/R-Run20221018_%j.out
#SBATCH --error=/work/lotterhos/MVP-NonClinalAF/slurm_log/R-Run20221018_%j.err

source ~/miniconda3/bin/activate MVP_env_R4.0.3

set -e
set -u
set -o pipefail

#### User modified values ####
# for 0b-final_params-fastruns-20220428.txt --array=2-1000%30
# for 0b-final_params-fastruns-20220428-b.txt #SBATCH --array=2-952%30
# for 0b-final_params-longruns-20220428.txt #SBATCH --array=2-301%30 #should be 301

# Local working path (this should navigate to the MVP repo)
mypath="/work/lotterhos/MVP-NonClinalAF"
cd ${mypath}

# Folder within MVP where you want are your output files
outpath="sim_output_20220428/"
runID="20220428"

# Parameter file
params="src/0b-final_params-fastruns-20220428-b.txt"
## Notes:
## For first set of fast runs, params="src/0b-final_params-fastruns-20220428.txt"
## For second set of fast runs, params="src/0b-final_params-fastruns-20220428-b.txt"
## For long runs, params="src/0b-final_params-longruns-20220428.txt"


# Extracting individual variables
# If you update or change your parameters see
# https://github.com/ModelValidationProgram/MVP-NonClinalAF/blob/alan/notebook/20210630_creatingSLiMBatchScriptVariables.md
# for guide on replacing this block of code.
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
    

##############
#### run R script  (takes a few  min)
#############
SECONDS=0 # used to time analyses
echo "Running R scripts"
Rscript --vanilla src/c-AnalyzeSimOutput.R ${seed} ${outpath} ${runID} > ${outpath}${seed}"_R.out" 2> ${outpath}${seed}"_R.error"
echo "Done with processing first R script. Analysis took $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min"


mv ${seed}"_genotypes.pcaProject" ${outpath}
