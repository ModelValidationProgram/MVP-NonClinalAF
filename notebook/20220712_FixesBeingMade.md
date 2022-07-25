# More minor fixes and bugs

### Slim Code

I found a minor bug in the simulation code, which will only affect the 1 trait simulations.
The normalization of fitness and fitness itself are calculated from  a difference S.D. of the fitness function:
```
 if (N_traits == 1)
		{
			fitness_norm = dnorm(0.0, 0.0, fitness_var_1);
			subpop.setValue("fitness_norm", fitness_norm);
		}
		fits = dnorm(asFloat(inds.getValue("phenotype1")), asFloat(opts[1]), fitness_var_2)/fitness_norm;
```

They should both be `fitness_var_2`

To avoid confusion, I renamed them to `fitness_var_0` and `fitness_var_1` and `SIGMA_K_0` and `SIGMA_K_1` so they correspond to `opt0` and `opt1` 
I made the changes and saved the Slim Code with the same filename.


### R Code

* `final.df$cor_af_sal_noutliers` and `temp`: it would be better to output the number of QTN outliers and the number of truly neutral outliers - I think this calculation migth be wrong
	*  The calculation is correct for `cor_af_temp_noutliers` and `cor_af_sal_noutliers`
	*  `cor_af_temp_noutliers = num_causal_sig_temp_corr + num_notCausal_sig_temp_corr` and same for salinity
		* the "noncausal" ones include neutral loci linked to selected loci
	* Here was the error:

WRONG CODE. This code was wrong because it was asking for NOT neutral loci, which would have been all loci on the first 10 LGs
```
num_neut_sig_temp_corr <- sum(muts_full$af_cor_temp_P[!muts_full$causal_temp=="neutral"]<Bonf_alpha)# number of neutral (unlinked to causal) loci that have significant  Kendall's correlations with temperature after Bonferroni correction
  num_neut_sig_sal_corr <- sum(muts_full$af_cor_sal_P[!muts_full$causal_sal=="neutral"]<Bonf_alpha)# number of neutral (unlinked to causal) loci that have significant Kendall's correlations with salinity after Bonferroni correction
```

CORRECTED CODE fixed 2022-07-21
```
num_neut_sig_temp_corr <- sum(muts_full$af_cor_temp_P[muts_full$causal_temp=="neutral"]<Bonf_alpha)# number of neutral (unlinked to causal) loci that have significant  Kendall's correlations with temperature after Bonferroni correction
  num_neut_sig_sal_corr <- sum(muts_full$af_cor_sal_P[muts_full$causal_sal=="neutral"]<Bonf_alpha)# number of neutral (unlinked to causal) loci that have significant Kendall's correlations with salinity after Bonferroni correction
```
  	  

* `LEA3.2_lfmm2_num_causal_sig_temp` there is a second variable `LEA3.2_lfmm2_num_causal_sig_temp` and analogous variable for salinity seems to be missing `LEA3.2_lfmm2_num_neut_sig_sal`
	* I fixed this 2022-07-21 

* add a check that the muts table matches the G file output
