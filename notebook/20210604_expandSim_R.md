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

- [ ] make sure can simulate confounding demography
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

### Plan to change R code (NEXT WEEK)
- [ ] Add info about the sim to the top of each plot
- [ ] Add mutation effect size to RDA plot

### Plan to test Sims  (NEXT WEEK)
- [ ] test some of the sims outlined below
- [x] test an obvious case that should produce clines
    - [ ] try reducing mutation rate and increasing effect size - see if can get bigger clines


### sims to test

`N_traits = 2, ispleiotropy=0, iscontrol=0, demog=Estuary` should give either salinity or temp mutations, but not mutations with effects on both. 
   -  [ ] check that no pleiotropy works - mutations have effects on one trait or the other
   
Decrease mutation rate and increase effect size for one cline and see if we get bigger effects


-  [ ]  `N_traits = 1, ispleiotropy=1, iscontrol=0, demog=SS` 
    *  I coded this so the simulation still simulates 2 traits, but only one of them (the latitude temp trait) is under selection. I haven't checked the common garden outputs to see if they make sense.
    *  I don't think it should matter what `ispleieotropy` is set to. If it is 1, mutations affect both traits, but only one of them is under selection. If it is 0, mutations affect one trait or the other, but only one trait is under selection.
    * this appears to be working. With `ispleieotropy=1`, I get an almost perfect correlation with the temperature trait and environment, and no correlation with the salinity trait and environment. 
    * this is run `2242330863068`
        * this sim wasn't perfect, I had high migration on x and low migration on y. But it should still produce strong clines.
    * So, it appears there are some temperature outliers, but also a large number of causal mutations. I might get better clines with a lower mutation rate
    - [ ] CHECK THAT THE OUTLIERS EXPLAIN A MAJOR PROPORTION OF VA
    

-  [ ]  `N_traits = 2, ispleiotropy=1, iscontrol=0, demog=Estuary` this is the situation I've been testing



-  [ ]  `N_traits = 2, ispleiotropy=0, iscontrol=0, demog=Estuary` should give either salinity or temp mutations, but not mutations with effects on both. 
    -  [ ] check that no pleiotropy works - mutations have effects on one trait or the other
    * For the non-pleiotropy case, I multiplied the QTN mutation rate by 2, because in the non-pleiotropy case I've hacked it so that the mutations are still pleiotropic (so the rest of the code works), but then I make them non-pleiotropic by randomly setting the mutational effect on one trait to 0. In essence, this halves the mutation rate for the non-pleiotropic case compared to the pleiotropic case for each trait. This hack (QTL_mu*2) make the mutation rates equal between the pleiotropic and non-pleiotropic cases.    


* then, test all the above with `iscontrol=1`

* Then, test all the above with `demog=SS` and equal migration on x and y


### Plan to set up parameter space
- [ ] set up a data frame of scenarios
- [ ] Tweak unix pipeline to run a scenario quickly
