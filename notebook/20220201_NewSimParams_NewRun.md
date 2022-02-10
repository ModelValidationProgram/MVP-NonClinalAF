
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
Submitted batch job 23041426
```

```
ls -l *_muts.txt | wc -l # this is the number of sims that completed the SliM run
999
ls -l *_plusneut_MAF01.recode2.vcf.gz | wc -l 
999

(base) [lotterhos@login-00 src]$ sbatch e-run_nonAF_sims_1R-fastruns-20220201.sh
Submitted batch job 23053326
```

#### Let's check Fst
awk 'NR==1' 1231094_Rout_simSummary.txt > ../summary_20220201.txt # header

awk 'FNR>1' *_Rout_simSummary.txt >> ../summary_20220201.txt # data

```
                                                   Fst        LA
SS-Clines_N-equal_m-constant                 0.03531793 0.4712589
SS-Mtn_N-equal_m-constant                    0.03889448 0.4654960
Est-Clines_N-equal_m-constant                0.07327978 0.4289701
SS-Clines_N-cline-N-to-S_m-constant          0.08639464 0.4660734
SS-Mtn_N-cline-N-to-S_m-constant             0.09169559 0.4598002
Est-Clines_N-cline-N-to-S_m-constant         0.13836340 0.4202269
SS-Clines_N-cline-center-to-edge_m-constant  0.08632038 0.4659329
SS-Mtn_N-cline-center-to-edge_m-constant     0.09422532 0.4608467
Est-Clines_N-cline-center-to-edge_m-constant 0.15824315 0.4229619
SS-Clines_N-equal_m_breaks                   0.04873994 0.4772067
SS-Mtn_N-equal_m_breaks                      0.05263310 0.4719693
Est-Clines_N-equal_m_breaks                         NaN       NaN
SS-Clines_N-variable_m-variable              0.07102395 0.4363785
SS-Mtn_N-variable_m-variable                 0.07894541 0.4109113
Est-Clines_N-variable_m-variable                    NaN       NaN
```

Old simulation parameters for SS in `0a-setUpSims.Rmd`
```
MIG_x = c(0.05, 0.05, 0.49)
MIG_y = c(0.05, 0.05, 0.1)
```

The first thing to notice is that I can not have similar levels of divergence and local adaptation in the SS vs Estuary demographies.

New simulation parameters:
```
MIG_x = c(0.03, 0.03, 0.49)
MIG_y = c(0.03, 0.03, 0.07)
```

This should bring FST of the SS closer to 0.05.

Old simulation parameters for N-variable m-variable:
```
  var_m_ss = c(rep(0.01,5),rep(0.1,10),rep(0.25,5));
```

New simulation params for N-variable m-variable:
```
  var_m_ss = c(0.001, rep(0.01,5),rep(0.1,10),rep(0.25,5));
```

### Failed sims

There are at least 65 failed sims.

``
1231193_genotypes.pcaProject  1231388_genotypes.pcaProject  1231603_genotypes.pcaProject  1231813_genotypes.pcaProject  1232003_genotypes.pcaProject  1232243_genotypes.pcaProject
1231198_genotypes.pcaProject  1231403_genotypes.pcaProject  1231613_genotypes.pcaProject  1231823_genotypes.pcaProject  1232008_genotypes.pcaProject
``


R. out file:
```
"LFMM finished"
<RDA output>
[1] 0.606522
[1] 0.7042434
[1] 0.1921053
[1] 0.8648148
List of 4
 $ x   : num [1:2] 0.7 1.9
 $ y   : num [1:2] 33.3 10.8
 $ xlab: NULL
 $ ylab: NULL
[1] 33.33877 10.83860
pdf
  2
pdf
  2
```

R.error file:
```
Warning message:
Removed 10 rows containing missing values (geom_point).
Error in sample.int(length(x), size, replace, prob) :
  cannot take a sample larger than the population when 'replace = FALSE'
Calls: sample -> sample.int
Execution halted
```

I'm willing to bet that some sims result in less than 3000 loci, and that's a problem.


Yes, I was correct in thinking that the number of loci was quite low in the N-variable m-variable scenario.

This led to errors in the RDA prediction, which needed 5000 loci.

![image](https://user-images.githubusercontent.com/6870125/152754989-0fc555f5-c1df-4b92-8f44-84bd7ea7089c.png)

## To do: Figure out how to "up" genetic diversity

Try a pyslim recaptiation with a larger number of individuals! Right now the "N" for pyslim is 1000, bring it up to 1000.

I'm going to be sloppy and overwrite the very first simulation.

Number of mutations produced by first run:
```
(base) [lotterhos@login-00 sim_output_20220201]$ wc -l 1231094_plusneut_MAF01.recode2.vcf.gz
45197 1231094_plusneut_MAF01.recode2.vcf.gz ## 10,000 individuals

(base) [lotterhos@login-00 sim_output_20220201]$ wc -l 1231094_Rout_Gmat_sample.txtTRUETRUE
15351 1231094_Rout_Gmat_sample.txtTRUETRUE # the 1,000 subsample individuals
```

#### N=1000 run
```
(base) [lotterhos@login-00 src]$ sbatch d-run_nonAF_sims_0Slim-fastruns-20220201-pyslimtest.sh
Submitted batch job 23075173
(base) [lotterhos@login-00 src]$ seff 23075173
```
This was still running after 7 hours.

#### N=200 run
```
(base) [lotterhos@login-00 src]$ sbatch d-run_nonAF_sims_0Slim-fastruns-20220201-pyslimtest.sh
Submitted batch job 23076368
```

(base) [lotterhos@login-00 src]$ seff 23076368
Job ID: 23076368
Array Job ID: 23076368_2
Cluster: discovery
User/Group: lotterhos/users
State: COMPLETED (exit code 0)
Cores: 1
CPU Utilized: 03:07:16
CPU Efficiency: 99.96% of 03:07:21 core-walltime
Job Wall-clock time: 03:07:21
Memory Utilized: 848.70 MB
Memory Efficiency: 16.58% of 5.00 GB


## The Rerun

- [x] rerun starting with R code to generate parameter list, Slim Code, and R analysis code

In summary, these simulations should bring FST of the SS scenarios to be higher, 
genetic diversity of all simulations to be higher (increased pyslim N for recap), 
and FST/diversity for the N-variable/m-variable simulations to be higher (added some lower migration rates)

```
(base) [lotterhos@login-00 src]$ sbatch d-run_nonAF_sims_0Slim-fastruns-20220201.sh
Submitted batch job 23086467
```

### Compare number of loci in vcf file for old vs. new

Old sims with N=100 recap: `45197 1231094_plusneut_MAF01.recode2.vcf.gz`

New sims with N=200 recap: `70193 1231094_plusneut_MAF01.recode2.vcf.gz`

Most will be filtered again after sampling.

The newer simulations take a lot longer to run. After 12 hours about 30-40% of 1000 runs are finished.

### Run R
```
(base) [lotterhos@login-01 src]$ sbatch e-run_nonAF_sims_1R-fastruns-20220201.sh
Submitted batch job 23128766
```

```
wc -l 1231094_Rout_Gmat_sample.txt
26092 1231094_Rout_Gmat_sample.txt
```

This is one of the N-variable m-variable simulations (still low diversity - probably due to too high gene flow):
```
(base) [lotterhos@login-01 sim_output_20220201]$ wc -l 1231103_Rout_Gmat_sample.txt
5538 1231103_Rout_Gmat_sample.txt
```

### Inspecting this round of sims
```
awk 'NR==1' 1231094_Rout_simSummary.txt > ../summary_20220201.txt # header
awk 'FNR>1' *_Rout_simSummary.txt >> ../summary_20220201.txt # data
```

```

                                                    Fst        LA
SS-Clines_N-equal_m-constant                 0.05782970 0.4841922
SS-Mtn_N-equal_m-constant                    0.06339836 0.4862318
Est-Clines_N-equal_m-constant                0.08777946 0.4297248
SS-Clines_N-cline-N-to-S_m-constant          0.12456721 0.4798241
SS-Mtn_N-cline-N-to-S_m-constant             0.13284979 0.4779571
Est-Clines_N-cline-N-to-S_m-constant         0.15853110 0.4226462
SS-Clines_N-cline-center-to-edge_m-constant  0.12918562 0.4781890
SS-Mtn_N-cline-center-to-edge_m-constant     0.13925447 0.4797316
Est-Clines_N-cline-center-to-edge_m-constant 0.18727872 0.4210846
SS-Clines_N-equal_m_breaks                   0.06927258 0.4886862
SS-Mtn_N-equal_m_breaks                      0.07520083 0.4905438
Est-Clines_N-equal_m_breaks                         NaN       NaN
SS-Clines_N-variable_m-variable              0.09337713 0.4365570
SS-Mtn_N-variable_m-variable                 0.10082130 0.4239514
Est-Clines_N-variable_m-variable                    NaN       NaN
```

![image](https://user-images.githubusercontent.com/6870125/153347014-f074d45e-1f75-4e34-93bf-870cc46394a1.png)
