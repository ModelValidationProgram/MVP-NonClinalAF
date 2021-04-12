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

-  [ ] check that no pleiotropy works - mutations have effects on one trait or the other
    * For the non-pleiotropy case, I multiplied the QTN mutation rate by 2, because in the non-pleiotropy case I've hacked it so that the mutations are still pleiotropic (so the rest of the code works), but then I make them non-pleiotropic by randomly setting the mutational effect on one trait to 0. In essence, this halves the mutation rate for the non-pleiotropic case compared to the pleiotropic case for each trait. This hack (QTL_mu*2) make the mutation rates equal between the pleiotropic and non-pleiotropic cases.    

- [ ] Situations to test:

* `N_traits = 1, ispleiotropy=1, iscontrol=0, demog=Estuary` (when N_traits = 1, set ispleieotropy to 1 because the way mutations are simulated. They will still have effects on the other phenotype, but that phenotype is not under selection)
*  `N_traits = 2, ispleiotropy=1, iscontrol=0, demog=Estuary` this is the situation I've been testing
*  `N_traits = 2, ispleiotropy=0, iscontrol=0, demog=Estuary` should give either salinity or temp mutations, but not mutations with effects on both. 

* then, test all the above with `iscontrol=1`
* Then, test all the above with `demog=SS` and equal migration on x and y


