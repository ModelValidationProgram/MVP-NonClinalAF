
# Revise Simulation parameters to give more realistic FST

### Edit `0a-setUpSims.Rmd`
Old simulation parameters:
```
c("SS-Clines", "SS-Mtn", "Est-Clines")
MIG_x = c(0.01, 0.01, 0.49)
MIG_y = c(0.01, 0.01, 0.01)
```

New simulation parameters in `0a-setUpSims.Rmd`:
```
MIG_x = c(0.1, 0.1, 0.49)
MIG_y = c(0.1, 0.1, 0.05)
```

(Need to run `0a-setUpSims.Rmd` with SimID 20220201)

### New simulation parameters in `a-PleiotropyDemog_20220101.slim`:

The m-breaks scenario and variable-m scenario:

Old code:
```
	if (demog=="Estuary"){

	var_m_estuary = c(0.001,rep(0.01,10),rep(0.1,10),rep(0.25,5));
	mig_breaks_estuary = 0.001;
```

New code:
```
	var_m_estuary = c(rep(0.01,5),rep(0.1,10),rep(0.25,5));
	mig_breaks_estuary = 0.01;
```

Old code:
```
	if (demog=="SS"){ // STEPPING STONE DEMOGRAPHY

	var_m_ss = c(rep(0.001,6),0.01,0.1,0.25);
  # mig_breaks was 0.0001
```  

New code:
```
  var_m_SS = c(0.005, rep(0.01,5),rep(0.1,10),rep(0.25,5));
	mig_breaks_SS = 0.01;
```

### Updating pipeline files

Updated files for a 20220201 simulation ID

```
(base) [lotterhos@login-01 src]$ sbatch d-run_nonAF_sims_0Slim-fastruns-20220201.sh
Submitted batch job 23002489
```

All sims finished:
```
ls -l *_muts.txt | wc -l # this is the number of sims that completed the SliM run
999 # all ran!

ls -l *_plusneut_MAF01.recode2.vcf.gz | wc -l 
999
```

Submit R code:
```
(base) [lotterhos@login-00 src]$ sbatch e-run_nonAF_sims_1R-fastruns-20220201.sh
Submitted batch job 23015482
```

Check runs finished:
```
ls *_pdf_1pop.pdf | wc -l # number of sims that were analyzed through the population step
999

ls -l *_2muts.pdf | wc -l # number that analyzed the mutations
999

ls -l *.pcaProject | wc -l # this is the number of sims that were analyzed through the PCA step in the R script
998

ls -l *_Rout_simSummary.txt | wc -l # this is the number of sims that were analyzed through the final output in the R script
998

awk 'NR==1' 1231094_Rout_simSummary.txt > ../summary_20220202.txt # header

awk 'FNR>1' *_Rout_simSummary.txt >> ../summary_20220202.txt # data
```

> TO DO: One simulation didn't finish. Figure out which one that is and fix it.

### Summary results

Fst and LA for each demography:
```
                                                   Fst        LA
SS-Clines_N-equal_m-constant                 0.01712748 0.4465919
SS-Mtn_N-equal_m-constant                    0.01895314 0.4286733
Est-Clines_N-equal_m-constant                0.10689793 0.4317783
SS-Clines_N-cline-N-to-S_m-constant          0.04712509 0.4394426
SS-Mtn_N-cline-N-to-S_m-constant             0.05045080 0.4236419
Est-Clines_N-cline-N-to-S_m-constant         0.18356182 0.4223818
SS-Clines_N-cline-center-to-edge_m-constant  0.04702098 0.4408296
SS-Mtn_N-cline-center-to-edge_m-constant     0.05146663 0.4251662
Est-Clines_N-cline-center-to-edge_m-constant 0.21886227 0.4259855
SS-Clines_N-equal_m_breaks                   0.03150117 0.4570688
SS-Mtn_N-equal_m_breaks                      0.03387413 0.4392155
Est-Clines_N-equal_m_breaks                         NaN       NaN
SS-Clines_N-variable_m-variable              0.72529261 0.4822981
SS-Mtn_N-variable_m-variable                 0.74568698 0.4806462
Est-Clines_N-variable_m-variable                    NaN       NaN
```

The Est-clines scenario has higher Fst, but similar level of local adaptation (at least within a demography).

The N-variable, m-variable scenarios still have too high an Fst, so I'll have to rerun those.

![image](https://user-images.githubusercontent.com/6870125/152299024-cd904f78-89a2-47e1-ba04-dc94a24cd3bd.png)

I think I should aim for an FST of 0.05, like we did in Lotterhos and Whitlock.

Number of Loci
```
                        oligogenic mod.\npolygenic highly\npolygenic
Pre MAF filter            12.24691       668.07385         3042.2980
Post filter MAF > 0 .01    7.62037        64.31077          525.5186
```
![image](https://user-images.githubusercontent.com/6870125/152300368-344c256c-d8eb-4f61-a273-39f7e168e5b9.png)

Main result is the same

![image](https://user-images.githubusercontent.com/6870125/152301019-f939785a-1439-4ba4-bfdd-4f3ada3e4371.png)

> TO DO: The simulation that has a negative correlation is 1231973. figure out why that is.

![image](https://user-images.githubusercontent.com/6870125/152301332-10b5b347-5e2d-4656-944c-1467bff576a7.png)

![image](https://user-images.githubusercontent.com/6870125/152301485-3ad37b30-bbbb-417e-b30c-d2ab750f730b.png)


# Revise Simulation parameters to give more realistic FST

### Edit `0a-setUpSims.Rmd`
Old simulation parameters in `0a-setUpSims.Rmd`:
```
MIG_x = c(0.1, 0.1, 0.49)
MIG_y = c(0.1, 0.1, 0.05)
```

New simulation parameters:
```
MIG_x = c(0.07, 0.07, 0.49)
MIG_y = c(0.07, 0.07, 0.07)
```

(Need to run `0a-setUpSims.Rmd` with SimID 20220201)

### New simulation parameters in `a-PleiotropyDemog_20220101.slim`:

The m-breaks scenario and variable-m scenario:

Old code:
```
	var_m_estuary = c(rep(0.01,5),rep(0.1,10),rep(0.25,5));
	mig_breaks_estuary = 0.01;
```

```
	if (demog=="SS"){ // STEPPING STONE DEMOGRAPHY
	var_m_ss = c(rep(0.001,6),0.01,0.1,0.25);
	mig_breaks_SS = 0.01;
```
WHOOPS! This was not changed as I intended.

New code:
```
	var_m_estuary = c(rep(0.01,5),rep(0.1,10),rep(0.25,5));
	mig_breaks_estuary = 0.01;
```
```
	if (demog=="SS"){ // STEPPING STONE DEMOGRAPHY
  var_m_SS = c(rep(0.01,5),rep(0.1,10),rep(0.25,5));
	mig_breaks_SS = 0.01;
```


## Reogranizing

```
mv sim_output_20220201/ archived/sim_output_20220201a`
mv summary_20220202.txt archived/summary_20220202a.txt
(base) [lotterhos@login-00 src]$ cp a-PleiotropyDemog_20220201.slim archived/a-PleiotropyDemog_20220201a.slim
(base) [lotterhos@login-00 src]$ cp 0b-final_params-fastruns-20220201.txt archived/0b-final_params-fastruns-20220201a.txt
```

I'm going to use the same simID, so I don't need to change all the code. I'm being lazy. The old outputs will be stored with the "a" at the end

```
(base) [lotterhos@login-00 src]$ sbatch d-run_nonAF_sims_0Slim-fastruns-20220201.sh
Submitted batch job 23023767
```

```

ls -l *_muts.txt | wc -l # this is the number of sims that completed the SliM run
846

ls -l *_plusneut_MAF01.recode2.vcf.gz | wc -l 
846

ls *_pdf_1pop.pdf | wc -l # number of sims that were analyzed through the population step
846

ls -l *_2muts.pdf | wc -l # number that analyzed the mutations
846

ls -l *.pcaProject | wc -l # this is the number of sims that were analyzed through the PCA step in the R script
846

ls -l *_Rout_simSummary.txt | wc -l # this is the number of sims that were analyzed through the final output in the R script
846

awk 'NR==1' 1231094_Rout_simSummary.txt > ../summary_20220202.txt # header

awk 'FNR>1' *_Rout_simSummary.txt >> ../summary_20220202.txt # data
```

ANOTHER BUG!!! This time in the SliM code. God help me. It's probably a semicolon.


### Summary results

Fst and LA for each demography:
```
                                             Fst        LA
SS-Clines_N-equal_m-constant                 0.02497782 0.4593779
SS-Mtn_N-equal_m-constant                    0.02760938 0.4492423
Est-Clines_N-equal_m-constant                0.08785785 0.4297248
SS-Clines_N-cline-N-to-S_m-constant          0.06532060 0.4543940
SS-Mtn_N-cline-N-to-S_m-constant             0.06931278 0.4454430
Est-Clines_N-cline-N-to-S_m-constant         0.15945164 0.4226462
SS-Clines_N-cline-center-to-edge_m-constant  0.06534542 0.4562008
SS-Mtn_N-cline-center-to-edge_m-constant     0.07142221 0.4454235
Est-Clines_N-cline-center-to-edge_m-constant 0.18708251 0.4210846
SS-Clines_N-equal_m_breaks                   0.03906303 0.4679972
SS-Mtn_N-equal_m_breaks                      0.04213524 0.4580187
Est-Clines_N-equal_m_breaks                         NaN       NaN
SS-Clines_N-variable_m-variable                     NaN       NaN
SS-Mtn_N-variable_m-variable                        NaN       NaN
Est-Clines_N-variable_m-variable                    NaN       NaN
```

FST still not high enough for the Baseline Clines, and too high for the Est Clines in the N-equal m-equal scenario

Old simulation parameters:
```
MIG_x = c(0.07, 0.07, 0.49)
MIG_y = c(0.07, 0.07, 0.07)
```

New simulation parameters:
```
MIG_x = c(0.05, 0.05, 0.49)
MIG_y = c(0.05, 0.05, 0.1)
```

Why did the N-equal m-breaks fail?

failed seeds - checking slim output `vi 1231103_outfile.error.txt`
```
Error on script line 157, character 51 (inside runtime script block):

               destSubpop.setMigrationRates(sourceID, sample(var_m_ss, 1));
                                                             ^^^^^^^^
							    ```
							    
I FOUND THE BUG!! It was a capital letter.

## Reogranizing

```
mv sim_output_20220201/ archived/sim_output_20220201b
mv summary_20220202.txt archived/summary_20220202b.txt
cp src/0b-final_params-fastruns-20220201.txt archived/0b-final_params-fastruns-20220201b.txt
```

I'm going to use the same simID, so I don't need to change all the code. I'm being lazy. The old outputs will be stored with the "b" at the end

## Next simulation set - hopefully this one gives FST = 0.05 for baseline scenario

* Reran `0a-setUpSims.Rmd` with new migration rates
```
(base) [lotterhos@login-01 src]$ sbatch d-run_nonAF_sims_0Slim-fastruns-20220201.sh
Submitted batch job 23040196
```
