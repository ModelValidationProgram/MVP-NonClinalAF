cd ../sim_output_150_20210830/ #update folder
awk 'NR==1' 1231094_Rout_simSummary.txt > ../summary.txt
awk 'FNR>1' *_Rout_simSummary.txt >> ../summary.txt
