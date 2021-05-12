
FINAL OUTPUTS
all_corr_phen_sal
subsamp_corr_phen_temp
all_corr_phen_temp
subsamp_corr_phen_temp
num_causal_muts
n_ind

MATRICES TO WRITE
af_pop



COMMON GARDEN STUFF
```{r, echo=FALSE, eval=FALSE}
# Sanity check
# the fitness column in the indPhen_df gives
# the fitnesses in the home population (subpop)
# these should match the fitnesses in the same
# subpop in the common garden data frame
# Individual common garden data
indCG_df<- read.table(paste0(path,seed,"_fitnessmat_ind.txt"))
dim(indCG_df)
rownames(indCG_df) = c("popID", 1:(nrow(indCG_df)-1))
head(indCG_df[,1:20])
# the first row in this is popID
head(indPhen_df)
head(indCG_df[,1:10])
cbind(round(indPhen_df$fitness[1:10],6), 
      unlist(round(indCG_df[2,1:10],6)))

# hard coding
cbind(round(indPhen_df$fitness[1990:2000],6), 
      unlist(round(indCG_df[21,1990:2000],6)))

cor(indPhen_df$fitness,indPhen_df$sal_opt)
cor(indPhen_df$phen_sal,indPhen_df$sal_opt)

```


```{r filter VCF and reformat}
# VCF file from SliM
#vcf <- read.vcfR(paste0(path,seed,"_VCF_causal.vcf.gz"))
#head(vcf)
#head(vcf@fix, 50)
#dim(vcf@fix)
# example of how to find a specific mutation in the vcf file
#  muts_df[1,]
#  vcf@fix[grep(muts_df$mutID[1], vcf@fix[,"INFO"]),]


# This website has code I can recycle here
# https://github.com/TestTheTests/TTT_RecombinationGenomeScans/blob/master/src/b_Proc_Sims.R

head(muts_df)
dim(muts_df)
dim(vcf@gt)
head(vcf@gt[,1:5])

head(vcf@fix)

geno <- vcf@gt[,-1] 
position <- getPOS(vcf)

if (sum(duplicated(position)) != 0){
  print("This simulation has duplicated locus positions")
}

G <- matrix(NA, nrow = nrow(geno), ncol = ncol(geno))
G[geno %in% c("0/0", "0|0")] <- 0
G[geno  %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
G[geno %in% c("1/1", "1|1")] <- 2

a_freq <- rowSums(G)/(2*ncol(G))
#hist(a_freq)
```

```{r figure out mtuations to keep}
vcf_muts0 <- stringr::str_extract_all(vcf@fix[,"INFO"], "MID=(\\d*);", simplify=TRUE)
head(vcf_muts0)
vcf_muts1 <- sub("MID=", "", vcf_muts0)
vcf_muts <- as.numeric(sub(";", "", vcf_muts1))
head(vcf_muts)
keepmuts <- which(vcf_muts %in% muts_df$mutID)
head(keepmuts)
```

```{r check vcf files line up}
dim(vcf@gt[keepmuts,])
dim(vcf_full@gt[keepmuts_full,])

if(!identical(as.numeric(position[keepmuts]), position_full[keepmuts_full]+1)){
  print("Error: positions in causal vcf vs. casual mutations in full vcf not the same")
}

if(!identical(vcf@gt[keepmuts,2], vcf_full@gt[keepmuts_full,2])){
  print("Error: genotypes in causal vcf vs. casual mutations in full vcf not the same")
}


if (!identical(sort(vcf_muts[keepmuts]), sort(muts_df$mutID))){
  print("Error: mutations in VCF not lining up with muts_df")
  break
}

```


The following code is for the subset VCF file of just causal mutations, subset to sampled individuals
```{r subset causal mutations MAF}
vcf_muts_subset <- vcf_muts[keepmuts]
muts_df$mutID

head(vcf@fix[keepmuts,], 50)

keepinds <- which(indPhen_df$subset)

dim(G)

G_subset <- G[keepmuts,keepinds]
dim(G_subset)
# should have 1000 individuals

rownames(G_subset)=vcf_muts_subset
colnames(G_subset)=keepinds
G_subset[,1:10]

## This is my subset of individuals
length(keepinds)
subset_indPhen_df <- indPhen_df[keepinds,]

#cbind(vcf@gt[,2], G[,1])
#cbind(vcf@gt[keepmuts,2], G[keepmuts,1]) # looks good 
```


```{r sanity check individual phenotypes}

## ADD INFO FROM SLIM VCF FILE
vcf@fix[keepmuts,"INFO"]
head(muts_df)
dim(muts_df)
for (i in 1:nrow(muts_df)){
  ind = grep(muts_df$mutID[i], vcf@fix[,"INFO"])
  if (length(ind)==1){
    muts_df$INFO[i] <-  vcf@fix[ind,"INFO"]
    muts_df$position[i] <- position[ind]
    muts_df$af_slimVCF[i] <- a_freq[ind]
  }
}


## ADD INFO FROM PYSLIM VCF FILE
for (i in 1:nrow(muts_df)){
  ind = grep(muts_df$mutID[i], vcf_full@fix[,"ALT"])
  if (length(ind)==1){
    muts_df$position_pyslim[i] <- position_full[ind]
    muts_df$af_pyslim[i] <- a_freq_full[ind]
  }
}

head(muts_df)

plot(muts_df$af_slimVCF, muts_df$af_pyslim)
abline(0,1)
```


ridge regression genotype-environment temp
```{r}
lfmm_ridge_ENVtemp <- lfmm::lfmm_ridge(Y = scaled.genotype, X = env_temp, K = 4, lambda = 1e-4)
#The lfmm.ridge object contains estimates for the latent variables and for the effect sizes. Here, the estimates are used for computing calibrated significance values and for testing associations between the response matrix Y and the explanatory variable x. It can be done as follows:
lfmm_ridge_ENVtemp_test  <- lfmm::lfmm_test(Y = scaled.genotype, X = env_temp, lfmm = lfmm_ridge_ENVtemp, calibrate = "gif")

lfmm_ridge_ENVtemp_P <- lfmm_ridge_ENVtemp_test$calibrated.pvalue

print(c("lfmm.test.ridge$gif", lfmm_ridge_ENVtemp_test$gif))
#hist(p.values.ridge, col = "lightgreen", main="LFMM ridge")
#qval <- qvalue::qvalue(p.values)
#plot(qval)
#The plot suggests that setting fdr.level = 0.025 warrant few false positives.
#qval <- qvalue::qvalue(p.values, fdr.level = 0.005)
#candidates <- which(qval$significant)
#plot(training$position, -log10(p.values.ridge), cex = .5, pch = 19, col = "black", main="LFMM ridge", ylim=c(0, 60))
#plot_layers(y_head=55, y_arrows=c(10, 0))

LFMM_ridge_ENVtemp_log10p <- -log10(as.numeric(lfmm_ridge_ENVtemp_P))
#plot(final_df$pos, final_df$LFMM_ridge_0.0_ALL_log10p) 

plot(muts_full$pos_pyslim, LFMM_ridge_ENVtemp_log10p)
points(muts_full$pos_pyslim[muts_full$causal], LFMM_ridge_ENVtemp_log10p[muts_full$causal], col="blue", pch=20)
```


ridge regression genotype-phenotype temp
```{r lfmm ridge temp phen}

lfmm_ridge_PHENtemp <- lfmm::lfmm_ridge(Y =  scaled.genotype, X = phen_temp, K = 4, lambda = 1e-4)
#The lfmm.ridge object contains estimates for the latent variables and for the effect sizes. Here, the estimates are used for computing calibrated significance values and for testing associations between the response matrix Y and the explanatory variable x. It can be done as follows:
lfmm_ridge_PHENtemp_test  <- lfmm::lfmm_test(Y = scaled.genotype, X = phen_temp, lfmm = lfmm_ridge_PHENtemp, calibrate = "gif")

lfmm_ridge_PHENtemp_P <- lfmm_ridge_PHENtemp_test$calibrated.pvalue

print(c("GIF", lfmm_ridge_PHENtemp_test$gif))
#hist(p.values.ridge, col = "lightgreen", main="LFMM ridge")
#qval <- qvalue::qvalue(p.values)
#plot(qval)
#The plot suggests that setting fdr.level = 0.025 warrant few false positives.
#qval <- qvalue::qvalue(p.values, fdr.level = 0.005)
#candidates <- which(qval$significant)
#plot(training$position, -log10(p.values.ridge), cex = .5, pch = 19, col = "black", main="LFMM ridge", ylim=c(0, 60))
#plot_layers(y_head=55, y_arrows=c(10, 0))

LFMM_ridge_PHENtemp_log10p <- -log10(as.numeric(lfmm_ridge_PHENtemp_P))
#plot(final_df$pos, final_df$LFMM_ridge_0.0_ALL_log10p) 

par(mar=c(4,4,1,1))
plot(muts_full$pos_pyslim, LFMM_ridge_PHENtemp_log10p)
points(muts_full$pos_pyslim[muts_full$causal], LFMM_ridge_PHENtemp_log10p[muts_full$causal], col="blue", pch=20)
```

```{r GWAS}
gwas_mut_slope <- NULL
gwas_mut_p <- NULL
gwas_mut_corrected_slope <- NULL
gwas_mut_corrected_p <- NULL

# this loop is slow but correct
for (i in 1:nrow(G_full_subset)){
  gwas_mut_i <- lm(G_full_subset[i,]~subset_indPhen_df$phen_temp)
  gwas_mut_slope[i] <- gwas_mut_i$coefficients[2]
  gwas_mut_p[i] <- summary(gwas_mut_i)$coefficients[2,4]
  
  gwas_mut_corrected_i <- lm(G_full_subset[i,]~subset_indPhen_df$phen_temp + pcaG$loadings[1,] + pcaG$loadings[2,])
  gwas_mut_corrected_slope[i] <- gwas_mut_corrected_i$coefficients[2]
  gwas_mut_corrected_p[i] <-  summary(gwas_mut_corrected_i)$coefficients[2,4]
}

plot(muts_full$pos_pyslim, -log10(gwas_mut_p))
plot(muts_full$pos_pyslim, -log10(gwas_mut_corrected_p))
points(muts_full$pos_pyslim[muts_full$causal],
       -log10(gwas_mut_corrected_p[muts_full$causal]), pch=20, col="magenta")

plot(LFMM_ridge_PHENtemp_log10p, -log10(gwas_mut_corrected_p))
points(LFMM_ridge_PHENtemp_log10p[muts_full$causal], -log10(gwas_mut_corrected_p)[muts_full$causal], col="magenta", pch=20)

plot(LFMM_ridge_ENVtemp_log10p, LFMM_ridge_PHENtemp_log10p)
points(LFMM_ridge_ENVtemp_log10p[muts_full$causal], LFMM_ridge_PHENtemp_log10p[muts_full$causal], col="magenta", pch=20)
```



```{r lfmm ridge sal}

env_sal <- scale(as.matrix(subset_indPhen_df$sal_opt))
# ridge regression

lfmm_ridge_sal <- lfmm::lfmm_ridge(Y = scaled.genotype, X = env_sal, K = 4, lambda = 1e-4)
#The lfmm.ridge object contains estimates for the latent variables and for the effect sizes. Here, the estimates are used for computing calibrated significance values and for testing associations between the response matrix Y and the explanatory variable x. It can be done as follows:
lfmm_ridge_sal_test  <- lfmm::lfmm_test(Y = scaled.genotype, X = env_sal, lfmm = lfmm_ridge_sal, calibrate = "gif")

lfmm_ridge_sal_P <- lfmm_ridge_sal_test$calibrated.pvalue

print(c("lfmm.test.ridge$gif", lfmm_ridge_sal_test$gif))
#hist(p.values.ridge, col = "lightgreen", main="LFMM ridge")
#qval <- qvalue::qvalue(p.values)
#plot(qval)
#The plot suggests that setting fdr.level = 0.025 warrant few false positives.
#qval <- qvalue::qvalue(p.values, fdr.level = 0.005)
#candidates <- which(qval$significant)
#plot(training$position, -log10(p.values.ridge), cex = .5, pch = 19, col = "black", main="LFMM ridge", ylim=c(0, 60))
#plot_layers(y_head=55, y_arrows=c(10, 0))

LFMM_ridge_sal_log10p <- -log10(as.numeric(lfmm_ridge_sal_P))
#plot(final_df$pos, final_df$LFMM_ridge_0.0_ALL_log10p) 

plot(muts_full$pos_pyslim, LFMM_ridge_sal_log10p)
points(muts_full$pos_pyslim[muts_full$causal], LFMM_ridge_sal_log10p[muts_full$causal], col="blue", pch=20)

plot(muts_full$Va_sal_prop, LFMM_ridge_sal_log10p)
plot(muts_full$Va_temp_prop, LFMM_ridge_ENVtemp_log10p)
```


```{r}
subset_indPhen_df$rand_env <- dmvnorm(m,
                                      mean=c(0,0), sigma=matrix(c(1, 0.5, 0.5, 1), ncol=2))
# This is a saddle, I think 

subset_indPhen_df$rand_env2 <- dmvnorm(m,
                                       mean=c(1,1), sigma=matrix(c(1, 0, 0, 1), ncol=2))
# This is like disease - it increases with salinity and temperature



plot(subset_indPhen_df$sal_opt,  subset_indPhen_df$temp_opt)
plot(subset_indPhen_df$sal_opt, subset_indPhen_df$rand_env)
plot(subset_indPhen_df$temp_opt, subset_indPhen_df$rand_env)
plot(subset_indPhen_df$sal_opt, subset_indPhen_df$rand_env2)
plot(subset_indPhen_df$temp_opt, subset_indPhen_df$rand_env2)

rdaout <- rda(t(G_subset)~ subset_indPhen_df$sal_opt + subset_indPhen_df$temp_opt + subset_indPhen_df$rand_env + subset_indPhen_df$rand_env2)
summary(rdaout)
screeplot(rdaout)
plot(rdaout)
scores <- scores(rdaout, choices=1:4)
loci.sc <- scores$species
ind.sc <- scores$sites
head(loci.sc)
head(ind.sc)
str(loci.sc)

plot(loci.sc[,1], loci.sc[,2])

head(subset_indPhen_df)
subset_indPhen_df$RDA1 <- ind.sc[,1]
subset_indPhen_df$RDA2 <- ind.sc[,2]
subset_indPhen_df$RDA3 <- ind.sc[,3]
plot(subset_indPhen_df$temp_opt ~ subset_indPhen_df$RDA1)
plot(scale(subset_indPhen_df$phen_temp) ~ scale(subset_indPhen_df$RDA1))
abline(0,1, col="red")

plot(subset_indPhen_df$sal_opt ~ subset_indPhen_df$RDA2)
plot(subset_indPhen_df$phen_sal ~ subset_indPhen_df$RDA2)

p <- subset_indPhen_df$RDA2*0.53*0.41 + subset_indPhen_df$RDA3*0.84*0.29
# predicted salinity phenotype is the RDA loading * eigenvalue of axis * loading of salinity on that axis
plot(subset_indPhen_df$sal_opt ~ p)
plot(scale(subset_indPhen_df$phen_sal) ~ scale(p))
abline(0,1, col="red")
```


```{r}
alpha <- 0.05/nrow(muts_all)
sig_temp <- muts_all$af_cor_temp_P<alpha
sig_sal <- muts_all$af_cor_sal_P<alpha

# Percent of VA explained by clinal patterns (e.g. sig. correlations)
sum(muts_all$Va_temp_prop[sig_temp]) # sanity check
sum(muts_all$Va_sal_prop[sig_sal]) # sanity check

sum(sig_temp) # number significant 
sum(sig_sal) # number significant

sum(sig_temp/nrow(muts_all))
sum(sig_sal/nrow(muts_all))

muts_all$FalseNegCorType <- NA
muts_all$FalseNegCorType[!sig_temp & !sig_sal] <- "Both_FN"
muts_all$FalseNegCorType[sig_temp & !sig_sal] <- "Temp_TP_Sal_FN"
muts_all$FalseNegCorType[!sig_temp & sig_sal] <- "Temp_FN_Sal_TP"
muts_all$FalseNegCorType[sig_temp & sig_sal] <- "Both_TP"

table(muts_all$FalseNegCorType)/nrow(muts_all)

# some visualization
plot(abs(muts_all$Va_temp_prop), abs(muts_all$af_slope_temp))
plot(abs(muts_all$Va_temp_prop), abs(muts_all$af_cor_temp))

plot(abs(muts_all$Va_sal_prop), abs(muts_all$af_slope_sal))
plot(abs(muts_all$Va_sal_prop), abs(muts_all$af_cor_sal))

hist(vcf_muts_df$af_cor_temp)
hist(vcf_muts_df$af_cor_sal)

hist(abs(vcf_muts_df$af_cor_temp))
hist(abs(vcf_muts_df$af_cor_sal))

hist(abs(vcf_muts_df$af_slope_temp))
hist(abs(vcf_muts_df$af_slope_sal))
```
