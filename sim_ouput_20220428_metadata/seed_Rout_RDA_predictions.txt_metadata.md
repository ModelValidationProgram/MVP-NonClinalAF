### seed_Rout_RDA_predictions.txt (output after analysis of simulations)

correlation between (a linear prediction of the weighted RDA loadings) and (the individual's temp phenotype)

* `nloci` the number of loci randomly sampled from the genome
* `Random_cor_temppredict_tempphen` the correlation between the RDA temperature prediction for an individual based on a random set of loci (see equation below) and the temperature phenotype of an individual
    * `temp_pred <- ind.sc[,1]*eigenvals(rdaout_rand)[1]*summary(rdaout_rand)$biplot[2,1] + ind.sc[,2]*eigenvals(rdaout_rand)[2]*summary(rdaout_rand)$biplot[2,2]`
* `Random_cor_salpredict_salphen` the correlation between the RDA salinity prediction for an individual based on a random set of loci (see equation below) and the salinity phenotype of an individual
    * `sal_pred <- ind.sc[,1]*eigenvals(rdaout_rand)[1]*summary(rdaout_rand)$biplot[1,1] + ind.sc[,2]*eigenvals(rdaout_rand)[2]*summary(rdaout_rand)$biplot[1,2]`
* `Random_cor_RDAtemppredict_tempphen_structcorr` same as the 2nd column, but for the RDA model with structure correction
* `Random_cor_RDAsalpredict_salphen_structcorr` same as the 3rd column, but for the RDA model with structure correction