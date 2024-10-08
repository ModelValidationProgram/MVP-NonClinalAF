## Sync files to git
```
git status
git add <added files in commit stages>
git commit
git pull
git push
```


## rerunning simulations

* Hopefully all the bugs are caught!
* Keeping old set of simulations with this ID until the new ones are run and cross-checked:
`mv sim_output_20220428 sim_output_20220428a`

## Run a revised set of simultions with RunID 20220428

Using the same RunID because only minor bug fixes, and less chance of introducing new bugs by changing the run ID

```
cd /work/lotterhos/MVP-NonClinalAF/src
sbatch d-run_nonAF_sims_0Slim-fastruns-20220428.sh #Submitted batch job 29260794
squeue -p lotterhos
sbatch d-run_nonAF_sims_0Slim-fastruns-20220428-b.sh 
sbatch d-run_nonAF_sims_0Slim-longruns-20220428.sh #batch job  29366383
```

## Check slim runs are finished

```
cd sim_output_20220428

ls -l *outfile.txt | wc -l # 2250!

ls -l *_muts.txt | wc -l # 2250!

ls -l *_plusneut_MAF01.recode2.vcf.gz | wc -l #2250!
```

## Run R scripts

```
sbatch e-run_nonAF_sims_1R-fastruns-20220428.sh #Submitted batch job 29399132

sbatch e-run_nonAF_sims_1R-fastruns-20220428-b.sh #Submitted batch job 29460812

sbatch e-run_nonAF_sims_1R-longruns-20220428.sh #Submitted batch job 29484119
```

## Check R runs are finished
```
ls -la *_pdf_1pop.pdf | grep "Jul" | wc -l #2250

ls -la *_Rout_simSummary.txt | grep "Jul" | wc -l #2250
```

additional checks (I overwrote these files, so just want to check all the files are recent and there aren't any old ones left)
```
cd sim_output_20220428
ls -la | grep "1231094" | wc -l #number of files expected for each sim = 50 files
ls -la | grep "May" | wc -l #should be 0 files # but it's 2250
ls -la | grep "Jul" | wc -l # 110250 # should be 112500, the missing are the 2250 files
ls -la | head -n 100 # each seed has one file from May and it is 1233339_genotypes.pca
rm -rf *_genotypes.pca
```
Now everything is up to date

## Create summary

```
cd /work/lotterhos/MVP-NonClinalAF/sim_output_20220428

awk 'NR==1' 1231094_Rout_simSummary.txt > ../summary_20220428_20220726.txt # header

awk 'FNR>1' *_Rout_simSummary.txt >> ../summary_20220428_20220726.txt # data
```

```
awk 'NR==1' 1231094_LA.txt > ../summary_LAthroughtime_20220428_20220726.txt # header

awk 'FNR>1' *_LA.txt >> ../summary_LAthroughtime_20220428_20220726.txt # data
```

```
awk 'NR==1' 1231094_Rout_RDA_predictions.txt | awk '{print $0, "seed"}' > ../summary_20220428_20220726_RDApredictions.txt # header

for i in {1231094..1233343} # may need to loop through seeds 1231094 to 1233344
do
awk 'FNR>1' $i"_Rout_RDA_predictions.txt" | awk -v myi=" $i" '{print $0 myi}'  >> ../summary_20220428_20220726_RDApredictions.txt # data
done
```
## Run LA script
* Open R Studio in OOD
* Run the code in `c2-ContributionToLA.R`

## Make graphs for publication
