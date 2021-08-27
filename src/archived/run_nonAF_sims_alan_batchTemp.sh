#!/bin/bash

set -e
set -u
set -o pipefail

mypath="/home/adowneywall/Github/MVP-NonClinalAF"
cd ${mypath}

#outpath="sim_outputs_testAlan/"

for i in {1..10} 
do
	# Storing parameter list
	#VARNAME=`sed "${i}q;d" src/0b-final_params.txt`
	# Extracting individual variables

	level=`awk NR==${i} src/0b-final_params.txt | awk '{print $1}'`
	reps=`awk NR==${i} src/0b-final_params.txt | awk '{print $2}'`
	arch=`awk NR==${i} src/0b-final_params.txt | awk '{print $3}'`
	# demog_name=`awk '{print $4}'`
	# demog_level_sub=`awk '{print $5}'`
	# demog_level=`awk '{print $6}'`
	# MIG_x=`awk '{print $7}'`
	# MIG_y=`awk '{print $8}'`
	# xcline=`awk '{print $9}'`
	# ycline=`awk '{print $10}'`
	# demog=`awk '{print $11}'`
	# METAPOP_SIDE_x=`awk '{print $12}'`
	# METAPOP_SIDE_y=`awk '{print $13}'`
	# Nequal=`awk '{print $14}'`
	# isVariableM=`awk '{print $15}'`
	# MIG_breaks=`awk '{print $16}'`
	# arch_level_sub=`awk '{print $17}'`
	# arch_level=`awk '{print $18}'`
	# MU_base=`awk '{print $19}'`
	# MU_QTL_proportion=`awk '{print $20}'`
	# SIGMA_QTN_1=`awk '{print $21}'`
	# SIGMA_QTN_2=`awk '{print $22}'`
	# SIGMA_K_1=`awk '{print $23}'`
	# SIGMA_K_2=`awk '{print $24}'`
	# N_traits=`awk '{print $25}'`
	# ispleiotropy=`awk '{print $26}'`
	# seed=`awk '{print $27}'`


echo param 1: ${level} param 2: ${reps} param 3: ${arch} >> temp.txt

done
