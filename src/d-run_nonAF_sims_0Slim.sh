#!/bin/bash

#SBATCH --job-name=SlimRun20210909
#SBATCH --mail-user=k.lotterhos@northeastern.edu
#SBATCH --mail-type=FAIL
#SBATCH --partition=lotterhos
#SBATCH --mem=10G
#SBATCH --nodes=1
#SBATCH --array=141%1
#SBATCH --output=/work/lotterhos/MVP-NonClinalAF/slurm_log/SlimRun20210909_%j.out
#SBATCH --error=/work/lotterhos/MVP-NonClinalAF/slurm_log/SlimRun20210909_%j.err

source ~/miniconda3/bin/activate MVP_env
# This is a CONDA environment I created on my own personal CONDA folder using the environment found in src/env/MVP_env.yml

set -e
set -u
set -o pipefail

#### User modified values ####

# Local working path (this should navigate to the MVP repo)
mypath="/work/lotterhos/MVP-NonClinalAF"
cd ${mypath}
# Folder within MVP where you want are your output files
outpath="sim_output_150_20210909/"
mkdir -p ${outpath} # make outpath directory if it doesn't exist

# Parameter file
params="src/0b-final_params_20210830.txt"

#### User variables ####
# N for pyslim
POP=100
# Minimum allele freq. for vcftools
MAF=0.01

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
#### run slim sims (takes 10 min)
#############

# Initial timestamp
echo "Running SLiM..."
timestamp_initial=`date`
echo ${timestamp_initial}

#run slim in background
slim -d MY_SEED1=${seed} -d MY_RESULTS_PATH1=\"${outpath}\" -d MIG_x1=${MIG_x} -d MIG_y1=${MIG_y} -d demog1=${demog} -d xcline1=${xcline} -d ycline1=${ycline} -d METAPOP_SIDE_x1=${METAPOP_SIDE_x} -d METAPOP_SIDE_y1=${METAPOP_SIDE_y} -d Nequal1=${Nequal} -d isVariableM1=${isVariableM} -d MIG_breaks1=${MIG_breaks} -d MU_base1=${MU_base} -d MU_QTL_proportion1=${MU_QTL_proportion} -d SIGMA_QTN_1a=${SIGMA_QTN_1} -d SIGMA_QTN_2a=${SIGMA_QTN_2} -d SIGMA_K_1a=${SIGMA_K_1} -d SIGMA_K_2a=${SIGMA_K_2} -d N_traits1=${N_traits} -d ispleiotropy1=${ispleiotropy} src/a-PleiotropyDemog_202109.slim > $outpath${seed}_outfile.txt 2> ${outpath}${seed}_outfile.error.txt    

## Final timestamp
echo "Done with SLiM sims"
time_stampFinishSlim=`date`
echo ${time_stampFinishSlim}
HOURS=$(echo $(date -d "$time_stampFinishSlim" +%H) - $(date -d "$timestamp_initial" +%H) | bc)
MINS=$(echo $(date -d "$time_stampFinishSlim" +%M) - $(date -d "$timestamp_initial" +%M) | bc)
SECS=$(echo $(date -d "$time_stampFinishSlim" +%S) - $(date -d "$timestamp_initial" +%S) | bc)
echo "Slim run took: ${HOURS} hours, ${MINS} minutes, and ${SECS} seconds"

gzip -f ${outpath}${seed}"_VCF_causal.vcf"

##############
#### run python script to process tree sequence results
#### Takes overnight on my laptop ~24 hours
#############

# Initial timestamp
echo "Processing tree sequences in python..."
time_stampTreeSeq=`date`
echo ${time_stampTreeSeq}
    
python3 src/b-process_trees.py -s ${seed} -r 1e-06 -mu ${MU_base} -N ${POP} -o ${outpath} > ${outpath}${seed}_pytree.out.txt 2> ${outpath}${seed}_pytree.error.txt
# when the SLiM simulation started, and the "N" is the N for the whole metapopulation

#Final timestamp
echo "Done with processing tree sequences."
time_stampFinishTreeSeq=`date`
echo ${time_stampFinishTreeSeq}
HOURS=$(echo $(date -d "$time_stampFinishTreeSeq" +%H) - $(date -d "$time_stampTreeSeq" +%H) | bc)
MINS=$(echo $(date -d "$time_stampFinishTreeSeq" +%M) - $(date -d "$time_stampTreeSeq" +%M) | bc)
SECS=$(echo $(date -d "$time_stampFinishTreeSeq" +%S) - $(date -d "$time_stampTreeSeq" +%S) | bc)
echo "Tree sequence run took: ${HOURS} hours, ${MINS} minutes, and ${SECS} seconds"

##############
#### Compress vcf files and filter for MAFs
#############

echo "Generating VCF file..."
time_stampVCF=`date`
echo ${time_stampVCF}

vcftools --vcf ${outpath}${seed}"_PlusNeuts".vcf --maf ${MAF} --out ${outpath}${seed}"_plusneut_MAF01" --recode --recode-INFO-all

# Fix initial formating bug in vcf code
# the sed lines fixes a bug with pyslim output. See 20210324_pyslim.md for info     
sed 's/\.		/\.	0	/g' ${outpath}${seed}"_plusneut_MAF01.recode.vcf" > ${outpath}${seed}"_plusneut_MAF01.recode2.vcf" 
    
gzip -f ${outpath}${seed}"_plusneut_MAF01.recode2.vcf"
    
rm ${outpath}${seed}"_plusneut_MAF01.recode.vcf"
rm ${outpath}${seed}"_PlusNeuts".vcf

## Script completes
echo "Script Complete."
timestamp_end=`date`
echo ${timestamp_end}
HOURS=$(echo $(date -d "$timestamp_end" +%H) - $(date -d "$timestamp_initial" +%H) | bc)
MINS=$(echo $(date -d "$timestamp_end" +%M) - $(date -d "$timestamp_initial" +%M) | bc)
SECS=$(echo $(date -d "$timestamp_end" +%S) - $(date -d "$timestamp_initial" +%S) | bc)
echo "Script took: ${HOURS} hours, ${MINS} minutes, and ${SECS} seconds"

##############
#### run R script  (takes 3 min)
#############
#SECONDS=0 # used to time analyses
#echo "Running R scripts"
#Rscript --vanilla ../src/b_Proc_Sims.R ${i} $simType > ${i}"_Invers_R.out" 2> ${i}"_Invers_R.error" & echo $!
#wait #wait until the last background process is finished
#echo "Done with processing first R script. Analysis took $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min"
