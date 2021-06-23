#!/bin/bash

#SBATCH --job-name=SlimTest
#SBATCH --mem=50G
#SBATCH --time=0-10:00:00
#SBATCH --mail-user=downey-wall.a@northeastern.edu
#SBATCH --mail-type=FAIL
#SBATCH --partition=lotterhos
#SBATCH --nodes=1
#SBATCH --output=/work/lotterhos/MVP-NonClinalAF/slurm_log/SlimTest_%j.out
#SBATCH --error=/work/lotterhos/MVP-NonClinalAF/slurm_log/SlimTest_%j.err

set -e
set -u
set -o pipefail

mypath="/work/lotterhos/MVP-NonClinalAF"
cd $mypath

outpath="sim_outputs/"
i=1231094 #this is MY_SEED

    cat > sim_outputs/1231094__oliogenic_1-trait__Est-Clines_N-cline-center-to-edge_m-constant.txt
    
    ##############
    #### run slim sims (takes 10 min)
    #############
    echo "running SLiM"
    SECONDS=0 # used to time analyses, no spaces around "=" sign
	#run slim in background
    slim -d MY_SEED1=1231094 -d "MY_RESULTS_PATH1='sim_outputs/'" -d MIG_x1=0.01 -d MIG_y1=0.01 -d "demog1='SS'" -d "xcline1='V'" -d "ycline1='linear'" -d METAPOP_SIDE_x1=10 -d METAPOP_SIDE_y1=10 -d Nequal1=1 -d isVariableM1=0 -d MIG_breaks1=0 -d MU_base1=1e-07 -d MU_QTL_proportion1=0.01 -d SIGMA_QTN_1a=0.4 -d SIGMA_QTN_2a=0.4 -d SIGMA_K_1a=0.5 -d SIGMA_K_2a=0.5 -d N_traits1=1 -d ispleiotropy1=0 src/a-PleiotropyDemog.slim > sim_outputs/1231094_outfile.txt 2> sim_outputs/1231094_outfile.error.txt #  & echo $!
    
    
    # NOTE: the "&" runs in the background on my laptop, but may not be needed on cluster
		
   # wait
    echo "Done with SLiM sims. Analysis took $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min"


    gzip -f ${outpath}${i}"_VCF_causal.vcf"

    ##############
    #### run python script to process tree sequence results
    #### Takes overnight on my laptop ~24 hours
    #############
    cd results
    echo "processing tree sequences in python"
    SECONDS=0 # used to time analyses
    
    # mu is MU_base
  python3 src/b-process_trees.py -s 1231094 -r 1e-06 -mu 1e-06 -N 1000 > sim_outputs/1231094_pytree.out.txt 2> sim_outputs/1231094_pytree.error.txt
    #python3 src/b-process_trees.py -s ${i} -r 1e-06 -mu 1e-06 -N 1000 > ${outpath}${i}"_pytree.out.txt" 2> ${outpath}${i}"_pytree.error.txt" ##& echo $!
    # If I understand recaptitaion correctly, it applies to the point in time prior to
    # when the SLiM simulation started, and the "N" is the N for the whole metapopulation
    # I'll use N=10000 for all sims, since that is close to the total metapop size for all sims
    # NOTE: the "&" runs in the background on my laptop, but may not be needed on cluster

    #wait #wait until the last background process is finished

    echo "Done with processing tree sequences. Analysis took $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min"


    ##############
    #### compress vcf files and filter for MAFs
    #############
    
    vcftools --vcf ${outpath}${i}"_PlusNeuts".vcf --maf 0.01 --out ${outpath}${i}"_plusneut_MAF01" --recode --recode-INFO-all
    
    
    # the sed lines fixes a bug with pyslim output. See 20210324_pyslim.md for info
    
    #sed -i 's/[.]\t\t/[.]\t0\t/g' ${outpath}${i}"_plusneut_MAF01.recode.vcf" #maybe will work in Linux
    
    sed 's/\.		/\.	0	/g' ${outpath}${i}"_plusneut_MAF01.recode.vcf" > ${outpath}${i}"_plusneut_MAF01.recode2.vcf" #this works on my mac
    # not sure if it will work on linux
    
    gzip -f ${outpath}${i}"_plusneut_MAF01.recode2.vcf"
    
    rm ${outpath}${i}"_plusneut_MAF01.recode.vcf"
    rm ${outpath}${i}"_PlusNeuts".vcf
    
    

    ##############
    #### run R script  (takes 3 min)
    #############
    #SECONDS=0 # used to time analyses
    #echo "Running R scripts"
    #Rscript --vanilla ../src/b_Proc_Sims.R ${i} $simType > ${i}"_Invers_R.out" 2> ${i}"_Invers_R.error" & echo $!
    #wait #wait until the last background process is finished
    #echo "Done with processing first R script. Analysis took $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min"
