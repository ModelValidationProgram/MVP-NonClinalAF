#!/bin/bash
set -e
set -u
set -o pipefail

mypath="/Users/lotterhos/Documents/GitHub/MVP-NonClinalAF/src"
cd $mypath
nreps=210 #number of reps
start=10900  #10900 10950
finish=$(($start + $nreps-1))
echo -e $start $finish


for i in $(seq $start $finish)
do
    cd $mypath
    echo -e "\n\n"
    echo $i
	outfile=${i}"_out.txt"
	outpath="../sim_outputs/"
	outpathfile=${outpath}${outfile}
    #echo $outpathfile

    ##############
    #### run slim sims (takes 10 min)
    #############
    echo "running SLiM"
    SECONDS=0 # used to time analyses, no spaces around "=" sign
	#run slim in background
	slim -d seed=$i -d mypath="'${outpath}'" ../src_sims/INSERTFILENAME.slim > $outpathfile 2> $outpathfile".error" & echo $!
		# input my seed and path to outputput
    wait
    echo "Done with SLiM sims. Analysis took $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min"

    # make sure output from last analysis exists
    if [ ! -f "results/${i}.trees" ]; then
        echo "trees file for ${i} not found, skipping"
        continue # skip rest of code and go to next $i
    fi
    
    gzip -f ${outpath}${i}"_VCF_causal.vcf"

    ##############
    #### run python script to process tree sequence results  (takes 1 min)
    #############
    cd results
    echo "processing tree sequences in python"
    SECONDS=0 # used to time analyses
    python3 ../src/a_process_trees.py -s ${i} > ${outpath}${i}"_pytree.out" 2> ${outpath}${i}"_pytree.error" & echo $!

	TO DO ADD MORE PARAMETERS AS INPUT TO PROCESS TREES

    wait #wait until the last background process is finished

    echo "Done with processing tree sequences. Analysis took $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min"


    ##############
    #### compress vcf files and filter for MAFs
    #############
    
    vcftools --vcf ${outpath}${i}"_PlusNeuts".vcf --maf 0.01 --out ${outpath}${i}"_plusneut_MAF01" --recode --recode-INFO-all
    
    
    # the sed lines fixes a bug with pyslim output. See 20210324_pyslim.md for info
    
    #sed -i 's/[.]\t\t/[.]\t0\t/g' ${outpath}${i}"_plusneut_MAF01.recode.vcf" #maybe will work in Linux
    
    sed 's/\.		/\.	0	/g' ${outpath}${i}"_plusneut_MAF01.recode.vcf" > ${outpath}${i}"_plusneut_MAF01.recode2.vcf" #this works on my mac
    
    gzip -f ${outpath}${i}"_plusneut_MAF01.recode2.vcf"
    
    rm ${outpath}${i}"_plusneut_MAF01.recode.vcf"
    rm ${outpath}${i}"_PlusNeuts".vcf
    
    

    ##############
    #### run R script  (takes 3 min)
    #############
    SECONDS=0 # used to time analyses
    echo "Running R scripts"
    Rscript --vanilla ../src/b_Proc_Sims.R ${i} $simType > ${i}"_Invers_R.out" 2> ${i}"_Invers_R.error" & echo $!
    wait #wait until the last background process is finished
    echo "Done with processing first R script. Analysis took $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min"
