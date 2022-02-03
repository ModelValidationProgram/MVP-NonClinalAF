
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

