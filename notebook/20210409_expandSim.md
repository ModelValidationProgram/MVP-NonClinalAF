# 20210409

### Plan: add options
- [x] Today I wanted to add some new options to the simulation:

```
defineConstant("N_traits", 2);   // number of traits (1 or 2)
defineConstant("ispleiotropy", 1); // a value of 0 for no pleiotropy, 1 for pleitropy
defineConstant("iscontrol", 0); // a value of 0 for stabilizing selection
defineConstant("demog", "Estuary");
        // "Estuary" for estuary demography
        // "SS" for stepping stone demography
defineConstant("Nequal", 3); // 0 for equal N, 1 for N cline on x-axis, 2 for N cline on y-axis, 3 for variable N
defineConstant("isVariableM", 1); // 0 for equal m, 1 for variable m
```


I also added a table `<seed>_popInfo.txt` for the information about each population, including the environment and xy location and population size

### Plan: to do
- [x] To do: add variable m and variable N
    * Done! It all checks out when plotted in R

- [x] check that the Ne cline, environment, and x-location makes sense
    * I wrote some R code to plot the demography from each sim. I'm happy with the results! Everything checks out. I also added code to simulate Ne variation in y-axis or x-axis

### Situations to test:

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


