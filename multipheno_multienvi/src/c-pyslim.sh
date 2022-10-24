#!/bin/bash

#SBATCH --job-name=SlimMultvar20221024
#SBATCH --mail-user=k.lotterhos@northeastern.edu
#SBATCH --mail-type=FAIL
#SBATCH --partition=lotterhos
#SBATCH --mem=10G
#SBATCH --nodes=1
#SBATCH --output=/work/lotterhos/MVP-NonClinalAF/slurm_log/SlimMultvar20221024.out
#SBATCH --error=/work/lotterhos/MVP-NonClinalAF/slurm_log/SlimMultvar20221024.err

source ~/miniconda3/bin/activate MVP_env
# This is a CONDA environment I created on my own personal CONDA folder using the environment found in src/env/MVP_env.yml

set -e
set -u
set -o pipefail

# Local working path (this should navigate to the MVP repo)
mypath="/work/lotterhos/MVP-NonClinalAF"
cd ${mypath}
# Folder within MVP where you want are your output files
outpath="multipheno_multienvi/output_multisim/"

##############
#### TREE SEQUENCING
#############

echo "Processing tree sequences in python..."
time_stampTreeSeq=`date`
echo ${time_stampTreeSeq}

python3 multipheno_multienvi/src/c-pyslim.py > ${outpath}_pytree.out.txt 2> ${outpath}_pytree.error.txt

# less ${outpath}_pytree.error.txt

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

vcftools --vcf ${outpath}"892657863_PlusNeuts".vcf --maf ${MAF} --out ${outpath}${seed}"_plusneut_MAF01" --recode --recode-INFO-all

# Fix initial formating bug in vcf code
# the sed lines fixes a bug with pyslim output. See 20210324_pyslim.md for info     
sed 's/\.		/\.	0	/g' ${outpath}${seed}"_plusneut_MAF01.recode.vcf" > ${outpath}${seed}"_plusneut_MAF01.recode2.vcf" 
    
gzip -f ${outpath}${seed}"_plusneut_MAF01.recode2.vcf"
    
rm ${outpath}${seed}"_plusneut_MAF01.recode.vcf"
rm ${outpath}${seed}"_PlusNeuts".vcf
