
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



```{r predict phen from RDA LFMM outliers}
sum(muts_full$LEA3.2_lfmm2_mlog10P_salenv_sig | muts_full$LEA3.2_lfmm2_mlog10P_tempenv_sig)

# temperature or salinity significant outliers
G_lfmm_out <- G_full_subset[muts_full$LEA3.2_lfmm2_mlog10P_salenv_sig | muts_full$LEA3.2_lfmm2_mlog10P_tempenv_sig,]

str(G_lfmm_out)

rdaout_lfmm <- rda(t(G_lfmm_out)~ subset_indPhen_df$sal_opt + subset_indPhen_df$temp_opt)
scores <- scores(rdaout_lfmm, choices=1:4)
loci.sc <- scores$species
ind.sc <- scores$sites
summary(rdaout_lfmm)$biplot

subset_indPhen_df$RDA_LFMMloci_temp_pred <- ind.sc[,1]*eigenvals(rdaout_lfmm)[1]*summary(rdaout_lfmm)$biplot[2,1] + ind.sc[,2]*eigenvals(rdaout_lfmm)[2]*summary(rdaout_lfmm)$biplot[2,2]

plot(scale(subset_indPhen_df$phen_temp), scale(subset_indPhen_df$RDA_LFMMloci_temp_pred))
abline(0,1)
# predicted temperature phenotype is the RDA loading * eigenvalue of axis * loading of temperature on that axis

(RDA_LFMMloci_cor_temppredict_tempphen <- cor(scale(subset_indPhen_df$phen_temp), scale(subset_indPhen_df$RDA_LFMMloci_temp_pred), method = "spearman"))

## Salinity ###
subset_indPhen_df$RDA_LFMMloci_sal_pred <- ind.sc[,1]*eigenvals(rdaout_lfmm)[1]*summary(rdaout_lfmm)$biplot[1,1] + ind.sc[,2]*eigenvals(rdaout_lfmm)[2]*summary(rdaout_lfmm)$biplot[1,2]

plot(scale(subset_indPhen_df$phen_sal), scale(subset_indPhen_df$RDA_LFMMloci_sal_pred))
abline(0,1)
# predicted salinity phenotype is the RDA loading * eigenvalue of axis * loading of salinity on that axis
(RDA_LFMMloci_cor_salpredict_salphen <-cor(scale(subset_indPhen_df$phen_sal), scale(subset_indPhen_df$RDA_LFMMloci_sal_pred), method = "spearman"))
```

```{r predict phen from RDA with RDA outliers}

sum(muts_full$RDA_mlog10P_sig)

# temperature or salinity significant outliers
G_rda_out <- G_full_subset[which(muts_full$RDA_mlog10P_sig),]
str(G_rda_out)

rdaout_rda <- rda(t(G_rda_out)~ subset_indPhen_df$sal_opt + subset_indPhen_df$temp_opt)
scores <- scores(rdaout_rda, choices=1:4)
loci.sc <- scores$species
ind.sc <- scores$sites
summary(rdaout_rda)$biplot

subset_indPhen_df$RDA_RDAloci_temp_pred <- ind.sc[,1]*eigenvals(rdaout_rda)[1]*summary(rdaout_rda)$biplot[2,1] + ind.sc[,2]*eigenvals(rdaout_rda)[2]*summary(rdaout_rda)$biplot[2,2]

plot(scale(subset_indPhen_df$phen_temp), scale(subset_indPhen_df$RDA_RDAloci_temp_pred))
abline(0,1)
# predicted temperature phenotype is the RDA loading * eigenvalue of axis * loading of temperature on that axis

(RDA_RDAloci_cor_temppredict_tempphen <- cor(scale(subset_indPhen_df$phen_temp), scale(subset_indPhen_df$RDA_RDAloci_temp_pred), method = "spearman"))

## Salinity ###
subset_indPhen_df$RDA_RDAloci_sal_pred <- ind.sc[,1]*eigenvals(rdaout_rda)[1]*summary(rdaout_rda)$biplot[1,1] + ind.sc[,2]*eigenvals(rdaout_rda)[2]*summary(rdaout_rda)$biplot[1,2]

plot(scale(subset_indPhen_df$phen_sal), scale(subset_indPhen_df$RDA_LFMMloci_sal_pred))
abline(0,1)
# predicted salinity phenotype is the RDA loading * eigenvalue of axis * loading of salinity on that axis
(RDA_RDAloci_cor_salpredict_salphen <-cor(scale(subset_indPhen_df$phen_sal), scale(subset_indPhen_df$RDA_RDAloci_sal_pred), method = "spearman"))
```


For multivariate trait prediction, we may want to predict how far the predicted trait value is from the true trait value. I'm not sure how interesting this is, given random loci are as good as using outliers in some cases.
```{r}

str(subset_indPhen_df)

for (i in 1:nrow(subset_indPhen_df)){
  # Distance between true phenotype and predicted from all loci
subset_indPhen_df$dist_phen_pred_allloci[i] <- 
  dist(rbind(
  cbind(scale(subset_indPhen_df$phen_sal)[i],   
        scale(subset_indPhen_df$phen_temp)[i]
        ), 
  cbind(scale(subset_indPhen_df$RDA_allloci_temp_pred)[i], 
        scale(subset_indPhen_df$RDA_allloci_sal_pred)[i]
  )
  )
)

# Distance between true phenotype and predicted from LFMM loci
subset_indPhen_df$dist_phen_pred_LFMMloci[i] <- dist(rbind(
  cbind(scale(subset_indPhen_df$phen_sal)[i],   
        scale(subset_indPhen_df$phen_temp)[i]
        ), 
  cbind(scale(subset_indPhen_df$RDA_LFMMloci_temp_pred)[i], 
        scale(subset_indPhen_df$RDA_LFMMloci_sal_pred)[i]
  )
  )
)

# Distance between true phenotype and predicted from RDA loci
subset_indPhen_df$dist_phen_pred_RDAloci[i] <-dist(rbind(
  cbind(scale(subset_indPhen_df$phen_sal)[i],   
        scale(subset_indPhen_df$phen_temp)[i]
        ), 
  cbind(scale(subset_indPhen_df$RDA_RDAloci_temp_pred)[i], 
        scale(subset_indPhen_df$RDA_RDAloci_sal_pred)[i]
  )
  )
)
}

par(mfrow=c(3,1))
hist(subset_indPhen_df$dist_phen_pred_allloci, xlim=c(0, 10))
hist(subset_indPhen_df$dist_phen_pred_LFMMloci, col=rgb(0,0,1,0.3), xlim=c(0,10))
hist(subset_indPhen_df$dist_phen_pred_RDAloci, col=rgb(1,0,1,0.3), xlim=c(0,10))
```



```{r NOT USING}
### Visualize genotypes at different temps ####
new_temp=subset_indPhen_df$temp_opt

#northern pops
heatmap(t(G_subset[a$colInd,new_temp==1]), Rowv = NA,  Colv = NA,
        main="High temp (north) genotypes",
        labRow = subset_indPhen_df$subpop[new_temp==1],
        xlab="Mutation ID",  cexCol = 0.3, 
        ylab="Low Sal<---Population--->High Sal")

heatmap(t(G_subset[a$colInd,new_temp==-0.111111]), Rowv = NA,  Colv = NA,
        main="Mid-latitude genotypes",
        labRow = subset_indPhen_df$subpop[new_temp==-0.111111],
        xlab="Mutation ID",  cexCol = 0.3,  
        ylab="Low Sal<---Population--->High Sal")

heatmap(t(G_subset[a$colInd,new_temp==-1]), Rowv = NA,  Colv = NA,
        main="Low temp (south) genotypes",
        labRow = subset_indPhen_df$subpop[new_temp==-1],
        xlab="Mutation ID",  cexCol = 0.3, 
        ylab="Low Sal<---Population--->High Sal")

```


```{r}

X <- t(G_subset)
Y <- cbind(indPhen_df$sal_opt[which(indPhen_df$subset)], indPhen_df$temp_opt[which(indPhen_df$subset)])
str(X)
str(Y)


set.seed(12380923)
rm <- sort(sample(indPhen_df$indID[-keepinds], 200, replace=FALSE)) # 200 random individuals
head(cbind(indPhen_df$indID[rm+1], rm))
rm <- rm+1 # because individuals start at 0, this can be an index now
indPhen_df$test <- FALSE
indPhen_df$test[rm] <- TRUE

indPhen_df[which(indPhen_df$subset & indPhen_df$test),]
  # shows no overlap

# unit test - make sure test and train individuals do not overlap
if (sum((indPhen_df$subset & indPhen_df$test))){
  print("Error: overlap in subset and test set")}

test_indPhen_df <- indPhen_df[rm,]

trainX <- X
trainY <- Y
testX <- t(G[keepmuts, rm])
testY_true <- cbind(indPhen_df$sal_opt[rm], indPhen_df$temp_opt[rm])
dim(testY_true)
testXsampName <- test_indPhen_df$indID

```


NEED TO EDIT THIS
```{r show cool patterns of fitness, fig.height=6, fig.width=4}
head(indPhen_df)
env_df <- unique(indPhen_df[,c("subpop","sal_opt","temp_opt")])
head(env_df)
(npops <- nrow(indCG_df)-1)

# choose an individual who has high fitness at cold end of range
# and in low salinity sites
head(indCG_df[,1:20], 10)
tail(indCG_df[,1:20], 10)

env_df$ind1_fit <- indCG_df[2:nrow(indCG_df), 1]
env_df$ind1_pch <- c(0,22)

# choose an individual who has high fitness at warm end of range
# and in low salinity sites
tail(indCG_df[,(ncol(indCG_df)-10):ncol(indCG_df)], 10)
env_df$ind2_fit <- indCG_df$V2000[2:nrow(indCG_df)]
env_df$ind2_pch <- c(1,21)

# choose an individual who has high fitness at cold end of range
# and in HIGH salinity sites
head(indCG_df[,101:120], 10)
env_df$ind3_fit <- indCG_df$V107[2:nrow(indCG_df)]
env_df$ind3_pch <- c(0,22)

# choose an individual who has high fitness at warm end of range
# and in HIGH salinity sites
tail(indCG_df[,(ncol(indCG_df)-10):ncol(indCG_df)], 10)
env_df$ind4_fit <- indCG_df$V1994[2:nrow(indCG_df)]
env_df$ind4_pch <- c(1,21)

### Start plot ####
par(mar=c(2,4,3,0.1), mfrow=c(2,1), oma=c(4,0,0,0))
### Plot low salinity adapted ####
  plot(env_df$temp_opt,env_df$ind1_fit, pch=env_df$ind1_pch, 
       col="blue", bg="grey", bty="l", las=1, ylim=c(0,1),
       xlab="", ylab="fitness", main="Low salinity adapted")
    points(ind1_fit~temp_opt, type="l", data=env_df, col=adjustcolor("blue",0.4))
    text(-1,1.08,"Ind. A", col="blue", adj=0)
  
  points(env_df$temp_opt,env_df$ind2_fit, pch=env_df$ind2_pch, col="red", bg="grey") 
  points(env_df$temp_opt,env_df$ind2_fit, pch=env_df$ind2_pch, col=adjustcolor("red", 0.2), type="l")
  text(0.6,1.08,"Ind. B", col="red", adj=0)
  
### Plot high salinity adapted ####
  plot(env_df$temp_opt,env_df$ind3_fit, pch=env_df$ind3_pch, 
       col="blue", bg="grey", bty="l", las=1, ylim=c(0,1),
       xlab="", ylab="fitness", main="High salinity adapted")
    points(ind3_fit~temp_opt, type="l", data=env_df, col=adjustcolor("blue",0.4))
    text(-1,1.08,"Ind. C", col="darkblue", adj=0)
  
  points(env_df$temp_opt,env_df$ind4_fit, pch=env_df$ind4_pch, col="red", bg="grey") 
  points(env_df$temp_opt,env_df$ind4_fit, pch=env_df$ind4_pch, col=adjustcolor("red", 0.2), type="l")
  text(0.6,1.08,"Ind. D", col="darkred", adj=0)

### Plot legend
par(xpd=NA)  
legend(-1,-0.4,  legend=c("high sal.", "low sal."), bty="n", adj=0, fill=c("grey", NA), horiz=TRUE)
mtext("temperature", side=1, outer=TRUE, adj=0.6)
```

```{r predict phen from RDA all loci}
subset_indPhen_df$RDA_allloci_temp_pred <- subset_indPhen_df$RDA1*eigenvals(rdaout)[1]*summary(rdaout)$biplot[2,1] + subset_indPhen_df$RDA2*eigenvals(rdaout)[2]*summary(rdaout)$biplot[2,2]

plot(scale(subset_indPhen_df$phen_temp), scale(subset_indPhen_df$RDA_allloci_temp_pred))
abline(0,1)
  # predicted temperature phenotype is the RDA loading * eigenvalue of axis * loading of temperature on that axis

(RDAallloci_cor_temppredict_tempphen <- cor(scale(subset_indPhen_df$phen_temp), scale(subset_indPhen_df$RDA_allloci_temp_pred), method = "spearman"))

## Salinity ###
subset_indPhen_df$RDA_allloci_sal_pred <- subset_indPhen_df$RDA1*eigenvals(rdaout)[1]*summary(rdaout)$biplot[1,1] + subset_indPhen_df$RDA2*eigenvals(rdaout)[2]*summary(rdaout)$biplot[1,2]

plot(scale(subset_indPhen_df$phen_sal), scale(subset_indPhen_df$RDA_allloci_sal_pred))
abline(0,1)
  # predicted salinity phenotype is the RDA loading * eigenvalue of axis * loading of salinity on that axis
(RDAallloci_cor_salpredict_salphen <-cor(scale(subset_indPhen_df$phen_sal), scale(subset_indPhen_df$RDA_allloci_sal_pred), method = "spearman"))
```

muts_full$af_cor_sal_pop = NA

colnames(af_sal) <- muts_full$mutID
rownames(af_sal) <- sal_levels
colnames(af_temp) <- muts_full$mutID
rownames(af_temp) <- temp_levels

for (row in 1:nrow(G_full_subset)){
    counts_sal <- table(G_full_subset[row,], subset_indPhen_df$sal_opt)
    counts_temp <- table(G_full_subset[row,], subset_indPhen_df$temp_opt)
      # this give a table of the number of alleles for individuals at that salinity
    forSum_sal <- counts_sal*as.numeric(rownames(counts_sal))
    forSum_temp <- counts_temp*as.numeric(rownames(counts_temp))
      # 0 for homozygote reference
      # 1 for heterozygote
      # 2 for homozygote derived
      # by multipling counts by number, heterozygotes are counted once
      # and homozygote derived are counted twice
      # also accounts for situations when only 2 genotypes found at a locus
  af_sal[,row] <- (colSums(forSum_sal))/(2*n_ind_sal)
  af_temp[,row] <- (colSums(forSum_temp))/(2*n_ind_temp)
  muts_full$af_cor_sal_pop[row] <- cor(af_sal[,row], sal_levels)
  muts_full$af_cor_temp_pop[row] <- cor(af_temp[,row], temp_levels)
  }


head(af_sal)