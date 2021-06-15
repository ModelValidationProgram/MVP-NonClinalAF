# 20210604


### Plan to change Sims (FRIDAY)
- [x] I changed the QTN mutation rate relative to the neutral mutation rate
    - it was 0.25*(neutral mu)
    - I changed it to something we can manipulate
    - e.g. 0.25 for polygenic, but maybe 0.01 for oligogenic
    - [ ] check if this still gives adaptation in the high migration scenario
- [x] add a neutral chromosome to the simulation
    - [x] I increasted the number of chromosomes from 10 to 20, but kept the QTNs on the first 10(?) chromosomes
    - [ ] check the number of chromosomes that are neutral
- [x] add some non-linear environments to the sims
    -  added `xcline` and `ycline` with options `linear` (-1 to 1) and `V` (-1 to 1 to -1)
    - [ ] check this gives good output
- [x] make sims more flexible for the two traits
    - if only one trait, SIGMA_K_1 is used
```
defineConstant("SIGMA_K_1", 0.5);                        // smaller is stronger stabilizing selection, // larger is weaker (wider) stabilizing selection for trait 1
defineConstant("SIGMA_K_2", 0.5);                        // smaller is stronger stabilizing selection, // larger is weaker (wider) stabilizing selection for trait 2
defineConstant("SIGMA_K_Cov", 0);                        // covariance in stabilizing selection
defineConstant("SIGMA_QTN_1", 0.1);            // standard deviation of mutational effect size - for trait 1
defineConstant("SIGMA_QTN_2", 0.1);            // standard deviation of mutational effect size - for trait 2
defineConstant("SIGMA_QTN_Cov", 0);            // covariance in mutation effect size (e.g. if negative, mutations that tend to have a positive effect on one trait, have a negative effect on the other)
```

- [x] add range center largest N
    - I added this scenario; it was hard to get the total metapopulation size to be exactly 10000, but I got it to be 10020

- [x] add biogeographic breaks
    - between populations 31-40 and 41- 50 
    - between populations 61-70 and 71-80
    - [ ] test if this works
```
defineConstant("MIG_breaks", 0) // 0 means no biogeographic breaks
    // a value of 1 introduces 2 biogeographic breaks
```

### Plan to set up parameter space
- [x] set up a data frame of scenarios
- [x ] Tweak unix pipeline to run a scenario quickly

### Issues
I fixed the code and started working on the command line arguments to get the sims to work.
```
slim-d MY_SEED1=1231094 -d "MY_RESULTS_PATH1='sim_outputs/'" -d MIG_x1=0.01 -d MIG_y1=0.01 -d "demog1='SS'" -d "xcline1='V'" -d "ycline1='linear'" -d METAPOP_SIDE_x1=10 -d METAPOP_SIDE_y1=10 -d Nequal1=4 -d isVariableM1=0 -d MIG_breaks1=0 -d MU_base1=1e-07 -d MU_QTL_proportion1=0.01 -d SIGMA_QTN_1a=0.4 -d SIGMA_QTN_2a=0.4 -d SIGMA_K_1a=0.5 -d SIGMA_K_2a=0.5 -d N_traits1=1 -d ispleiotropy1=0src/a-PleiotropyDemog.slim > sim_outputs/1231094_outfile.txt 2> sim_outputs/1231094_outfile.txt.error.txt
```
This parameter set will run for a few hundred gens (200-400 depending on the seed) and then crap out. I traced the error to the N_traits=1 fitness function.
I think I fixed it (some phenotypes were integers instead of numeric), but I need to double check which trait opt0 corresponds to. Right now the fitness function for Trait 1 corresponds to N_traits=1

Yup, the salinity trait is locally adapted, but I wanted it to be the temperature trait! I fixed it. 

### To do
- [x ] fix R code to write command line argument equal to what works in the notes above
- [x] edit SLiM code to make the temperature trait adapt for N_traits=1, and add more comments to code for opt0 and opt1
- [x] pyslim - 
I started running pyslim on my laptop for sim 1231094__oliogenic_1-trait__Est-Clines_N-cline-center-to-edge_m-constant.txt, with -N 10000 and -mu 1e-07, but it was taking forever (in 18 hours it wasn't finished). So I kept Ne*mu constant and changed it to  -N 1000 and -mu 1e-06, to see if it would run faster. That took about 24 hours on my laptop.


### Checks and sim tests
- [x] 1-trait sims seem to be working as expected (0 salinity effect, temp effect)
- [x] 2-trait sims with no pleitropy also seem to be working as expected (if a mutation has an effect, it is non-zero for 1 trait and 0 for the other trait)
- [x] oligogenic scenario - producing about 20 mutations - will have to look at R evaluation outputs, but I think probably only a few mutations are contributing the most to adaptation
    - "oliogenic_1-trait__Est-Clines_N-cline-center-to-edge_m-constant\n1231094"
- [x] the large N in the range center and low N at the range edge is working correctly
    - "oliogenic_1-trait__Est-Clines_N-cline-center-to-edge_m-constant\n1231094"
- [x] the SS-mtn gives an West-to-East environment cline
    - "1231453__oliogenic_2-trait-pleiotropy-equal-S__SS-Mtn_N-variable_m-variable"
    - this gave fairly high levels of LA for both traits (cor phen-env > 0.9)
    - the N-variable m-variable doesn't produce very interesting structure. I'm also worried drift doesn't really operate non-uniformly because of all the migration. I'm wondering if I should sample from smaller M's more often to create some breaks, so I changed the code. The new code still produces high levels of LA, just more breaks and sometimes isolated populations. After I run the prelim sims, I might sample from even smaller m's to create more breaks.
- [x] the SS-Est oligogenic 2 traits with pleiotropy gives adaptation to salinity after sampling
    - "oliogenic_2-trait-pleiotropy-equal-S__Est-Clines_N-equal_m-constant" , 1231142, 
    - this gives ~25 mutations with mu_proportion equal to 0.01 and sigma_qtn = 0.4. The cor_sal_ind was 0.81, so pretty high. I decreased the mu_proportion to be 0.001 to see what would happen. There was still high correlation with the salinity trait and salinity, which was great! This gave ~14 mutations. So I decreased by an order of magnitude for all simulations.
    
- [x] the SS-Est polygenic 2 traits with pleiotropy gives adaptation to salinity after sampling 
    - "polygenic_2-trait-pleiotropy-equal-S__Est-Clines_N-equal_m-constant" , 1231217
    - noticed an issue in the migration rates in the upper right hand corner of the sim - FIXED
    
- [x] the m-breaks works
    - "oliogenic_2-trait-pleiotropy-equal-S__Est-Clines_N-equal_m_breaks", 1231141
        - noticed an issue with the m_breaks not working properly - FIXED ( I think I understand how migration works in SLiM now!)
    - "oliogenic_2-trait-pleiotropy-equal-S__SS-Clines_N-equal_m_breaks", 1231146
        - looks great!
    
- [x] the unqual-S works
    - "oliogenic_2-trait-pleiotropy-unequal-S__Est-Clines_N-variable_m-variable", 1231158
        - there was definitely a lower correlation between the temp trait and temperature, I think it was the demog.
        - this was a cool sim because there was a lot of variation in "how locally adapted" each population was!
        - was run with `sigma_K2a=2`
    - "oliogenic_2-trait-pleiotropy-unequal-S__Est-Clines_N-equal_m-equal", 1231157
        - at first with `sigma_K2a=2`, this gave similar correlations between temp and trait then with equally strong S, so I changed `sigma_K2a=4`
        - this definitely gave less adaptation in the temperature trait, but still a good amount of LA


### Plan to change R code 
- [x] output table of sims and load it with R code
- [ ] add plot of cor(phen, env) through time
- [ ] incorporate R code into bash script
- [ ] Add info about the sim to the top of each plot
- [ ] Add mutation effect size to RDA plot

### Next steps
- [ ] send list of 150 prelim sims to Alan, he can benchmark running sims, pyslim, and vcftools on cluster
- [ ] Alan will write array to run the 150 sims
- [ ] after they run, I will spot check the outputs with R code
- [ ] run the R code on 150 sims
- [ ] summarize the results of the 150 sims
- [ ] scale up and run the whole she-bang
