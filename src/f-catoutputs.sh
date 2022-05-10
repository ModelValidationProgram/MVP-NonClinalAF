cd ../sim_output_20220428/ #update folder

ls -l *_muts.txt | wc -l # this is the number of sims that completed the SliM run

ls -l *_plusneut_MAF01.recode2.vcf.gz | wc -l # this is the number of sims that completed the vcf filtering

ls -l ../*.pcaProject | wc -l # this is the number of sims that were analyzed through the PCA step in the R script

ls -l *_Rout_simSummary.txt | wc -l # this is the number of sims that were analyzed through the final output in the R script



awk 'NR==1' 1231094_Rout_simSummary.txt > ../summary_20220117.txt # header

awk 'FNR>1' *_Rout_simSummary.txt >> ../summary_20220117.txt # data
