
LFMM trait prediction

https://bcm-uga.github.io/lfmm/articles/lfmm#prediction-of-phenotypic-values-with-polygenic-risk-scores 
```{r}

## Simulated phenotypes for Arabidopsis thaliana SNP data
data("example.data")
## Simulated (and real) methylation levels for sun-exposed tissue sampled
data("skin.exposure")
Y <- example.data$genotype
dim(Y)
X <- example.data$phenotype #scaled phenotype
dim(X)

mod.lfmm <- lfmm_ridge(Y = Y, 
                       X = X, 
                       K = 6)

pv <- lfmm_test(Y = Y, 
                X = X, 
                lfmm = mod.lfmm, 
                calibrate = "gif")

pvalues <- pv$calibrated.pvalue 

sig <- which(p.adjust(pvalues, method="fdr") < 0.25)
sig #<- which(pvalues< 0.01)


pred <- predict_lfmm(Y = Y, 
                     X = X,
                     fdr.level = 0.1, 
                     mod.lfmm)

length(pred$candidates)
plot(pred$prediction ~ X)

```