# 20210507

### Plan
- [x] finishing outputing outliers for LFMM and for RDA
- [x] test an obvious case that should produce clines
    - [ ] try reducing mutation rate and increasing effect size - see if can get bigger clines

- [ ] add a neutral chromosome to the simulation
- [ ] add some non-linear environments to the sims
- [ ] make sims more flexible for the two traits

- [ ] test some of the sims outlined below

- [ ] formalize all outputs
- [ ] add color to RDA plots to show the temp and salinity optimums of individuals - maybe use green to blue for salinity and light to dark for temperature
- [ ] fix heat maps so the colors are constrained

- [ ] set up a data frame of scenarios
- [ ] Tweak unix pipeline to run a scenario quickly


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


### current thoughts

Based on what I've looked at so far, it looks like LFMM has a high false negative rate but few false positives.
RDA has a high false positive rate.

However, we can put their outliers into an RDA and still make a pretty accurate prediction about the optimal environment for an individual.
So even though we can't identify the genetic basis of adaptation, in some cases we can still make an accurate prediction about the optimal environment for that genotype, 
and that prediction is robust to different sets of loci used.
