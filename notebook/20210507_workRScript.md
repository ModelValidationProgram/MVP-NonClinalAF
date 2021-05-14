# 20210507

### Plan for R code
- [x] finishing outputing outliers for LFMM and for RDA
- [x] formalize all outputs
- [x] add color to RDA plots to show the temp and salinity optimums of individuals - maybe use green to blue for salinity and light to dark for temperature

- [ ] fix heat maps so the colors are constrained
- [ ] finish all output for R code and summary stats for a sim


pdf(paste0(path,seed,"_pdf_muts.pdf"), width=5, height=5)

ggplot(pop_df) + geom_point(aes(x=x, y=y,size=opt0)) + geom_point(aes(x=x, y=y, color=opt1)) + theme_classic() + scale_colour_gradient2(high=rgb(1,0.4,0.2), low="cornflowerblue", mid=rgb(0.8,0.8,0.7), name="Temp") + geom_text(aes(x=x, y=y+0.3,label=N)) + labs(size="salinity") + ggtitle(paste0("Optimums; N traits = ", info$Ntraits))


ggplot(subset_indPhen_df) + geom_point(aes(x=PC1, y=PC2,size=sal_opt)) + geom_point(aes(x=PC1, y=PC2, color=temp_opt)) + theme_classic() + scale_colour_gradient2(high=rgb(1,0.4,0.2), low="cornflowerblue", mid=rgb(0.8,0.8,0.7), name="Temp") + labs(size="salinity") + ggtitle(paste0(seed,"; N traits = ", info$Ntraits))


### Plan to change Sims
- [ ] add a neutral chromosome to the simulation
- [ ] add some non-linear environments to the sims
- [ ] make sims more flexible for the two traits
- [ ] make sure can simulate confounding demography

### Plan to change R code

### Plan to test Sims
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


Paul Rawson
Damian Brady local bouys
New Meadows River
seeing if we can get approx at other locations
Most low Salinity 24-26 ppt
Not much freshwater
UM reports of oyster plants and experiments - begins to show complexity
Cascoe Bay - lots of oyster culture
Great Bay - Nature Conservancy Site
Wells might be a good place, wild
Sheepscott - monitoring site, DMR - we need to work with them - Marcy Nelson - 
Dana shellfish officer for the state - she might know if anyone collects from Cascoe 
Paul can see if people can help us in other places, but not Great Bay

### current thoughts

Based on what I've looked at so far, it looks like LFMM has a high false negative rate but few false positives.
RDA has a high false positive rate.

However, we can put their outliers into an RDA and still make a pretty accurate prediction about the optimal environment for an individual.
So even though we can't identify the genetic basis of adaptation, in some cases we can still make an accurate prediction about the optimal environment for that genotype, 
and that prediction is robust to different sets of loci used.
