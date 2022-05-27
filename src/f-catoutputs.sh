cd ../sim_output_20220428/ #update folder

ls -l *_muts.txt | wc -l # this is the number of sims that completed the SliM run

ls -l *_plusneut_MAF01.recode2.vcf.gz | wc -l # this is the number of sims that completed the vcf filtering

ls -l ../*.pcaProject | wc -l # this is the number of sims that were analyzed through the PCA step in the R script

ls -l *_Rout_simSummary.txt | wc -l # this is the number of sims that were analyzed through the final output in the R script


## Full summary

awk 'NR==1' 1231094_Rout_simSummary.txt > ../summary_20220428.txt # header

awk 'FNR>1' *_Rout_simSummary.txt >> ../summary_20220428.txt # data


## RDA subsets



awk 'NR==1' 1231094_Rout_RDA_predictions.txt | awk '{print $0, "seed"}' > ../summary_20220428_RDApredictions.txt # header
for i in {1231094..1233343} # may need to loop through seeds 1231094 to 1233344
do
awk 'FNR>1' $i"_Rout_RDA_predictions.txt" | awk -v myi=" $i" '{print $0 myi}'  >> ../summary_20220428_RDApredictions.txt # data
done	

#less ../summary_20220428_RDApredictions.txt


# get first 7 characters of filename (seed) tmp=${filename:0:7}

# add column 

