


### Compare to a regression with the environmental variable as the respone variable and PCs of the genotype matrix as the predictor variable

```{r}
library(LEA)
write.lfmm(G, "genotypes.lfmm")
pc = pca("genotypes.lfmm", 10, scale = TRUE)
str(pc)
plot(pc, xlim=c(1,10))
# Choose first 2 PC axes for model

PCmod <- lm(ind$env1_mat ~ pc$projections[,1] + pc$projections[,2] + pc$projections[,3])
summary(aov(PCmod))
plot(ind$phenotype1_mat, predict(PCmod))
plot(predict(PCmod), MATtraitPredict)
cor(predict(PCmod), ind$phenotype1_mat)
cor(MATtraitPredict, ind$phenotype1_mat)

## Trait 2

PCmod2 <- lm(ind$env2_MTWetQ ~ pc$projections[,1] + pc$projections[,2] + pc$projections[,3])
summary(aov(PCmod2))
cor(predict(PCmod2), ind$phenotype2_MTWetQ)
plot(predict(PCmod2), ind$phenotype2_MTWetQ)

cor(rda_trait_pred(rdaout, 2, K),  ind$phenotype2_MTWetQ)
plot(rda_trait_pred(rdaout, 2, K),  ind$phenotype2_MTWetQ)

## Trait 3
cor(rda_trait_pred(rdaout, 3, K),  ind$phenotype3_MTDQ)
plot(rda_trait_pred(rdaout, 3, K),  ind$phenotype3_MTDQ)
PCmod3 <- lm(ind$env3_MTDQ ~ pc$projections[,1] + pc$projections[,2] + pc$projections[,3])
summary(aov(PCmod3))
cor(predict(PCmod3), ind$phenotype3_MTDQ)
plot(predict(PCmod3), ind$phenotype3_MTDQ)

## Trait 4
cor(rda_trait_pred(rdaout, 4, K),  ind$phenotype4_PDM)
plot(rda_trait_pred(rdaout, 4, K),  ind$phenotype4_PDM)
PCmod4 <- lm(ind$env4_PDM ~ pc$projections[,1] + pc$projections[,2] + pc$projections[,3])
summary(aov(PCmod4))
cor(predict(PCmod4), ind$phenotype4_PDM)
plot(predict(PCmod4), ind$phenotype4_PDM)

## Trait 5
cor(rda_trait_pred(rdaout, 5, K),  ind$phenotype5_PwarmQ)
plot(rda_trait_pred(rdaout, 5, K),  ind$phenotype5_PwarmQ)
PCmod5 <- lm(ind$env5_PwarmQ~ pc$projections[,1] + pc$projections[,2] + pc$projections[,3])
summary(aov(PCmod5))
cor(predict(PCmod5), ind$phenotype5_PwarmQ)
plot(predict(PCmod5), ind$phenotype5_PwarmQ)

## Trait 6
cor(rda_trait_pred(rdaout, 5, K),  ind$phenotype5_PwarmQ)
plot(rda_trait_pred(rdaout, 5, K),  ind$phenotype5_PwarmQ)
PCmod5 <- lm(ind$env5_PwarmQ~ pc$projections[,1] + pc$projections[,2] + pc$projections[,3])
summary(aov(PCmod5))
cor(predict(PCmod5), ind$phenotype5_PwarmQ)
plot(predict(PCmod5), ind$phenotype5_PwarmQ)

```



