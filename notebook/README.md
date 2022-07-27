## Sync files to git
```
git status
git add <added files in commit stages>
git commit
git pull
git push
```


## Run a set of simulations with RunID 20220428

```
cd /work/lotterhos/MVP-NonClinalAF/src
sbatch d-run_nonAF_sims_0Slim-fastruns-20220428.sh
squeue -p lotterhos
sbatch d-run_nonAF_sims_0Slim-fastruns-20220428-b.sh
sbatch d-run_nonAF_sims_0Slim-longruns-20220428.sh
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
sbatch e-run_nonAF_sims_1R-fastruns-20220428.sh

sbatch e-run_nonAF_sims_1R-fastruns-20220428-b.sh

sbatch e-run_nonAF_sims_1R-longruns-20220428.sh
```

## Check R runs are finished
```
ls *_pdf_1pop.pdf | wc -l

ls -l *_Rout_simSummary.txt | wc -l
```

## Create summary
```
awk 'NR==1' 1231094_Rout_simSummary.txt > ../summary_20220428.txt # header

awk 'FNR>1' *_Rout_simSummary.txt >> ../summary_20220428.txt # data
```

## Run LA script
* Open R Studio in OOD
* Run the code in `c2-ContributionToLA.R`

## Make graphs for publication
* `g-FinalAnalysis.R`
