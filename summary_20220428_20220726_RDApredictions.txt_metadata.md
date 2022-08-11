## `summary_20220428_20220726_RDApredictions.txt` 

A table summarizing the accuracy of the RDA-predicted trait value.

-------------------------
### Field descriptions
-------------------------
* `nloci` Number of loci randomly sampled from the genome and then used to run the RDA
* `Random_cor_RDAtemppredict_tempphen` For this random sample of loci, the correlation between the RDA-predicted unscaled temperature trait value and the ground-truth evolved temperature phenotype. The RDA model did not correct for structure.
* `Random_cor_RDAsalpredict_salphen` For this random sample of loci, the correlation between the RDA-predicted unscaled Env2 trait value and the ground-truth evolved Env2 phenotype. The RDA model did not correct for structure.
* `Random_cor_RDAtemppredict_tempphen_structcorr` For this random sample of loci, the correlation between the pRDA-predicted unscaled temperature trait value and the ground-truth evolved temperature phenotype. The pRDA model included a structure correction.
* `Random_cor_RDAsalpredict_salphen_structcorr` For this random sample of loci, the correlation between the pRDA-predicted unscaled Env2 trait value and the ground-truth evolved Env2 phenotype. The pRDA model included a structure correction.
* `seed` simulation seed