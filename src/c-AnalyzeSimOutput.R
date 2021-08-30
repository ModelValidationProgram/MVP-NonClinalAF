# 
# 
# 
# packages_needed <- c("vcfR", "distances","ggplot2",  "fields", "stringr", "vegan", "robust", "mvtnorm", "viridis", "gridExtra", "devtools")
# 
# for (i in 1:length(packages_needed)){
#   if(!(packages_needed[i] %in% installed.packages())){install.packages(packages_needed[i])}
# }
# 
# if (!requireNamespace("BiocManager", quietly = TRUE)){install.packages("BiocManager")}
# 
# if (!requireNamespace("LEA", quietly = TRUE)){BiocManager::install("LEA")} 
# 
# if (!requireNamespace("qvalue", quietly = TRUE)){
#   BiocManager::install("qvalue")
# }
# 
# if(!requireNamespace("ggExtra", quietly=TRUE)){
#   devtools::install_github("daattali/ggExtra")
# }
# 
# 
# if(!requireNamespace("OutFLANK", quietly=TRUE)){
#   devtools::install_github("whitlock/OutFLANK")
# }



libraries_needed <- c("vcfR", "distances","ggplot2",  "fields", "stringr", "vegan", "robust", "mvtnorm", "viridis", "gridExtra", "devtools", 
                     "MLmetrics", "PRROC", "qvalue", "OutFLANK", "BiocManager", "LEA", "ggExtra")
for (i in 1:length(libraries_needed)){
  library( libraries_needed[i], character.only = TRUE, lib.loc = "/home/lotterhos/R/x86_64-pc-linux-gnu-library/4.0") # for OOD
#  library(libraries_needed[i], character.only = TRUE, lib.loc = "/home/lotterhos/miniconda3/envs/MVP_env_R4.0.3/lib/R/library") # for bash script
}


### load params ####

args = commandArgs(trailingOnly=TRUE)

#seed = 1231222
#path = "../sim_output_150_20210826/"

seed = args[1]
path = args[2]

### load data ####
vcf_full <- read.vcfR(paste0(path,seed,"_plusneut_MAF01.recode2.vcf.gz"))

  # Individual phenotype and fitness data in home pop
  indPhen_df <- read.table(paste0(path,seed,"_ind.txt"), header=TRUE, 
                           colClasses = c("character", rep("numeric",8)))
  
  pop_df <- read.table(paste0(path,seed,"_popInfo.txt"), header=TRUE, 
                       colClasses = c("character", rep("numeric",6))) 
  
  # Local adaptation through time df
  LA_df <- read.table(paste0(path,seed,"_LA.txt"), header=TRUE, 
                      colClasses = c("character", rep("numeric",10)),
                      na.strings = "NAN")
  

  muts_df <- read.table(paste0(path,seed,"_muts.txt"), header=TRUE, 
                        colClasses = c("character","numeric","character", rep("numeric",5)))
  
  vatemp <- muts_df$p*(1-muts_df$p)*muts_df$mutTempEffect^2
  vasal <- muts_df$p*(1-muts_df$p)*muts_df$mutSalEffect^2
  #plot(round(vatemp/sum(vatemp),2), abs(muts_df$cor_temp), xlim=c(0,1), ylim=c(0,1))
  vasal <- muts_df$p*(1-muts_df$p)*muts_df$mutSalEffect^2
  #if(sum(vasal)>0){
  #  plot(round(vasal/sum(vasal),2), abs(muts_df$cor_sal), xlim=c(0,1), ylim=c(0,1))
  #}
  
  info <- read.table(paste0(path,seed,"_info.txt"), header=TRUE, 
                     colClasses = c("character", 
                                    rep("numeric",15),
                                    "character",
                                    rep("numeric",6)))
  
  
  allsims <- read.table("0b-final_params.txt", header=TRUE)
  thissim <- allsims[grep(seed, allsims$seed),]
  (plotmain <- paste(thissim$level, seed, sep="\n"))

#### plot subpops and migration ####
  
  mig <- read.table(paste0(path,seed,"_popInfo_m.txt"), header=TRUE)
  
  par(mfrow=c(1,1), mar=c(4,4,3,1))
  pdf(paste0(path,seed,"_pdf_1pop.pdf"), width=6, height=6)
  
  if (var(mig$m)>0){
    mig_thick <- log10(mig$m) + 5
  }else{
    mig_thick <- rep(1, length(mig$m))
  }
  plot(pop_df$x, pop_df$y, col = rgb(0,0,0,0), bty="l", xlab="x", ylab="y", main=plotmain, cex.main=0.5)
  for (i in 1:nrow(mig)){
    start_x <- pop_df$x[which(pop_df$subpopID==mig$from[i])]
    start_y <- pop_df$y[which(pop_df$subpopID==mig$from[i])]
    end_x <- pop_df$x[which(pop_df$subpopID==mig$to[i])]
    end_y <- pop_df$y[which(pop_df$subpopID==mig$to[i])]
    adj = 0.01  
    if (end_x < start_x){end_x <- end_x +0.2; start_x <- start_x - 0.5}
    if (end_x > start_x){end_x <- end_x -0.2; start_x <- start_x + 0.5}
    if (end_y < start_y){end_y <- end_y +0.2; start_y <- start_y - 0.5}
    if (end_y > start_y){end_y <- end_y -0.2; start_y <- start_y + 0.5}
    arrows(start_x,start_y,end_x, end_y, col="cornflowerblue", lwd=mig_thick[i], length=0.1)
  }
  text(pop_df$x, pop_df$y, pop_df$N, cex=0.6)
  
  ggplot(pop_df) + geom_point(aes(x=x, y=y,size=opt0)) + geom_point(aes(x=x, y=y, color=opt1)) + theme_classic() + scale_colour_gradient2(high=rgb(1,0.4,0.2), low="cornflowerblue", mid=rgb(0.8,0.8,0.7), name="Temp") + geom_text(aes(x=x, y=y+0.3,label=subpopID)) + labs(size="Env2") + ggtitle(plotmain) + theme(plot.title = element_text(size = 7))
  
  plot(LA_df$gen, LA_df$cor_temp_ind, type="l", col="orange", lwd=3, main=plotmain, 
       ylab="cor (Env, phenotype)", xlab="Generation", ylim=c(0,1))
  
  if (sum(LA_df$cor_sal_ind, na.rm=TRUE)>0){
    points(LA_df$gen, LA_df$cor_sal_ind, col="cornflowerblue", type="l", lwd=3, lty=2)
  }
  legend(0,1, c("Temp", "Env2"), lwd=3, col=c("orange", "cornflowerblue"), lty=c(1,2), bty="n")
  
  plot(LA_df$gen, LA_df$local_adapt, type="l", col="cornflowerblue", main=plotmain, 
       ylab="Amount of local adaptation", xlab="Generation", ylim=c(0,1))
  
  plot(LA_df$gen, LA_df$mean_phen1, type="l", col="cornflowerblue", main=plotmain, 
       ylab="Mean metapopulation phenotype", xlab="Generation", ylim=c(-1,1))
    if (sum(LA_df$cor_sal_ind, na.rm=TRUE)>0){
      points(LA_df$gen, LA_df$mean_phen0, col="orange", type="l", lwd=3, lty=2)
    }
  
  # Find high fitness individuals ####
  npops <- length(levels(factor(indPhen_df$subpop)))
  n = 10 # number of individual per pop
  subset <- c()
  for (i in 1:npops){
    #set.seed(139982)
    bob <- sample(indPhen_df$indID[indPhen_df$subpop==i], size = n,
                  replace=FALSE, prob = indPhen_df$fitness[indPhen_df$subpop==i])
    subset <- c(subset, bob)
  }
  length(subset)
  head(subset)
  indPhen_df$subset <- FALSE
  indPhen_df$subset[sort((subset+1))] <- rep(TRUE, length(subset))
  # add one to subset here because indID starts at 0
  
  nind_samp <- sum(indPhen_df$subset)

  indPhen_df$subpopID <- indPhen_df$subpop
  indPhen_df <- merge(indPhen_df, pop_df)
  indPhen_df <- indPhen_df[order(indPhen_df$indID),]
  # check everyone in correct order:
  # plot(indPhen_df$indID)
  
  boxplot(indPhen_df$fitness~indPhen_df$subpop, ylab = "fitness", xlab="pop", main="all", ylim=c(0,1))
  boxplot(indPhen_df$fitness~indPhen_df$N, ylab = "fitness", xlab="N", main="all", ylim=c(0,1))
  
  boxplot(indPhen_df$fitness[indPhen_df$subset]~indPhen_df$subpop[indPhen_df$subset], ylab = "fitness", xlab="pop", main="subset", ylim=c(0,1))
  
  boxplot(indPhen_df$fitness[indPhen_df$subset]~indPhen_df$N[indPhen_df$subset], xlab="N", ylab="fitness")
  # looks good to me
  dev.off()
  ## end plot high fitness individuals ####
  
  if (var(indPhen_df$phen_sal)>0){
    all_corr_phen_sal <- cor(indPhen_df$phen_sal,indPhen_df$sal_opt)
    subsamp_corr_phen_sal <- cor(indPhen_df$phen_sal[which(indPhen_df$subset)], indPhen_df$sal_opt[which(indPhen_df$subset)])
  }else{
    all_corr_phen_sal <- NA
    subsamp_corr_phen_sal <- NA
  }
  
  all_corr_phen_temp <- cor(indPhen_df$phen_temp,indPhen_df$temp_opt)
  
  subsamp_corr_phen_temp <- cor(indPhen_df$phen_temp[which(indPhen_df$subset)],indPhen_df$temp_opt[which(indPhen_df$subset)])

### show cool patterns of mutation ####
  
  pdf(paste0(path,seed,"_pdf_2muts.pdf"), width=5, height=5)
  
  plot(muts_df$mutSalEffect, muts_df$mutTempEffect, ylim=c(-0.7, 0.7), xlim=c(-0.7,0.7), xlab="mutation effect on Env2 phenotype", ylab="mutation effect on temperature phenotype", bty="l", main=plotmain)
  abline(h=0, col="grey")
  abline(v=0, col="grey")
  
  geno_full <- vcf_full@gt[,-1] 
  dim(geno_full)
  position_full <- getPOS(vcf_full)
  rown <- as.numeric(vcf_full@fix[,"ALT"]) # mutation ID
  # rown = 1 is a neutral mutation
  causal_mut_locs <- which(rown %in% muts_df$mutID) # causal mutations
  #position_full[causal_mut_locs] # positions of causal mutations
  
  
  head(geno_full[,1:5])
  
  
  G_full <- matrix(NA, nrow = nrow(geno_full), ncol = ncol(geno_full))
  G_full[geno_full %in% c("0/0", "0|0")] <- 0
  G_full[geno_full  %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
  G_full[geno_full %in% c("1/1", "1|1")] <- 2
  rm(geno_full)
  
  
  #head(G_full[,1:5])
  
  a_freq_full <- rowSums(G_full)/(2*ncol(G_full))
  hist(a_freq_full, breaks=seq(0,1,0.01))
  
  if (!identical(sort(unique(unlist(as.numeric(G_full)))), as.numeric(0:2))){
    print("Error: full genotype matrix not uniquely 012")
    break
  }
  
  
  G_full_subset <- G_full[,which(indPhen_df$subset)]
    # really important here that the order of individuals in indPhen_df is in the same order as G_full_subset
   # In a line of code above, I made sure that the order was correct
  rm(G_full)
  rownames(G_full_subset)=vcf_full@fix[,"POS"]
  colnames(G_full_subset)=indPhen_df$indID[which(indPhen_df$subset)]
  head(G_full_subset[,1:5])
  plot(colnames(G_full_subset), indPhen_df$indID[which(indPhen_df$subset)])
  
  
  a_freq_subset <- as.numeric(rowSums(G_full_subset)/(2*ncol(G_full_subset)))
  # should have 1000 individuals
  #hist(a_freq_subset)
  
  muts_full0 <- data.frame(seed = as.character(seed),
                           VCFrow = 1:nrow(vcf_full@gt),
                           mutID = vcf_full@fix[,"ALT"],
                           pos_pyslim = as.integer(vcf_full@fix[,"POS"]),
                           a_freq_full = a_freq_full,
                           a_freq_subset = a_freq_subset)
  dim(muts_full0)
  
  muts_full <- merge(muts_full0, muts_df, by=c("mutID", "seed"), all.x=TRUE)
  muts_full <- muts_full[order(muts_full$pos_pyslim),]
  # make sure mutations in correct order
  plot(muts_full$pos_pyslim)
  rm(muts_full0)
  rm(muts_df)
  dim(muts_full)
  
  muts_full$mutTempEffect <- as.numeric(muts_full$mutTempEffect)
  muts_full$mutSalEffect <- as.numeric(muts_full$mutSalEffect)
  muts_full$causal <- FALSE
  muts_full$causal[muts_full$mutID!=1] <- TRUE
  
  
  # FILTER OUT ALLELS WITH MAF < 0.01 AFTER SAMPLING ####
   if(!identical(as.numeric(muts_full$pos_pyslim), as.numeric(rownames(G_full_subset)))){print("Error 1a: mutations not lined up");break()}
    
    rm_loci <- which(muts_full$a_freq_subset < 0.01 & !muts_full$causal)
  
    numCausalLowMAFsample <- sum(muts_full$a_freq_subset < 0.01 & muts_full$causal)
    
    # Filter out MAF loci from G_full_subset
    
    if(length(rm_loci)>0){
      G_full_subset <- G_full_subset[-rm_loci,]
      muts_full <- muts_full[-rm_loci,]
    }
    if(!identical(as.numeric(muts_full$pos_pyslim), as.numeric(rownames(G_full_subset)))){print("Error 1b: mutations not lined up");break()}
    
    dim(G_full_subset)
    head(G_full_subset[,1:10])
    
    hist(muts_full$a_freq_subset, breaks=seq(0,1,0.01))
  plot(muts_full$a_freq_full, muts_full$a_freq_subset)
  hist(muts_full$mutTempEffect, xlim=c(-1,1), xlab="Mutation temp effect", breaks=seq(-1,1,0.01), main=plotmain)
  hist(muts_full$mutSalEffect, xlim=c(-1,1), xlab="Mutation Env2 effect", breaks=seq(-1,1,0.01), main=plotmain)
  
 
  
  # Sanity check that muts_full mutations line up with G_full_subset
  if(nrow(muts_full)!=nrow(G_full_subset)){print("Error 1: mutation data frames have different number of rows");break()}
  
  head(muts_full$VCFrow)
  #Sanity check
  if(cor(muts_full$VCFrow, muts_full$pos_pyslim)<0.999){print("Error 2: order of mutations wrong");break()} #sanity check
  

  num_causal <- sum(muts_full$causal)
  num_neut <- sum(!muts_full$causal)
  
  dim(G_full_subset)
  
  subset_indPhen_df <- indPhen_df[indPhen_df$subset,]
  if(!identical(subset_indPhen_df$indID, sort(subset_indPhen_df$indID))){print("Error individuals in wrong order"); break()}
  rm(indPhen_df)
  str(subset_indPhen_df)
  
  head(pop_df)
  

  
  af_pop <- matrix(NA, nrow=length(subpop_subset), ncol=nrow(G_full_subset))

  
  n_ind <- table(subset_indPhen_df$subpop)
  
  subset_temp_opt <- tapply(subset_indPhen_df$temp_opt,
                            subset_indPhen_df$subpop, mean)
  
  subset_sal_opt <- tapply(subset_indPhen_df$sal_opt,
                           subset_indPhen_df$subpop, mean)
  
  subpop_subset <- as.numeric(rownames(subset_sal_opt))
  rownames(af_pop) <- subpop_subset
  colnames(af_pop) <- rownames(G_full_subset)
  
  muts_full$af_cor_temp <- NA
  muts_full$af_slope_temp <- NA
  muts_full$af_cor_temp_P <- NA
  muts_full$af_cor_sal <- NA
  muts_full$af_slope_sal <- NA
  muts_full$af_cor_sal_P <- NA
  

  # Calc AF on subset of individuals sampled ####
  
  for (row in 1:nrow(G_full_subset)){
    counts <- table(G_full_subset[row,], subset_indPhen_df$subpop)
    # this give a table of the number of alleles for individuals at that subpopulation
    
    forSum <- counts*as.numeric(rownames(counts))
    # 0 for homozygote reference
    # 1 for heterozygote
    # 2 for homozygote derived
    # by multipling counts by number, heterozygotes are counted once
    # and homozygote derived are counted twice
    # also accounts for situations when only 2 genotypes found at a locus
    af_pop[,row] <- (colSums(forSum))/(2*n_ind)
    
    # temp cors
    corrtemp <- cor.test(af_pop[,row], subset_temp_opt, method="kendall")
    muts_full$af_cor_temp[row] <- corrtemp$estimate
    muts_full$af_cor_temp_P[row] <- corrtemp$p.value
    muts_full$af_cor_temp_mlog10P[row] <- -log10(corrtemp$p.value)
    muts_full$af_slope_temp[row] <- lm(af_pop[,row]~subset_temp_opt)$coef[2]
    
    ## do it for salinity / Env2
    corrsal <- cor.test(af_pop[,row], subset_sal_opt, method="kendall")
    muts_full$af_cor_sal[row] <- corrsal$estimate
    muts_full$af_cor_sal_P[row] <- corrsal$p.value
    muts_full$af_cor_sal_mlog10P[row] <- -log10(corrsal$p.value)
    muts_full$af_slope_sal[row] <- lm(af_pop[,row]~subset_sal_opt)$coef[2]
    
  }
  
  # Sanity check that what I calculate here (from the subsample)
  # is similar to what I outputted from SliM in muts_df
  
  plot( muts_full$cor_temp, muts_full$af_cor_temp, ylim=c(-1,1), xlim=c(-1,1), xlab="SliM cor(p, Temp)", ylab="Subsample cor(p, Temp)", main=plotmain)
  abline(0,1)
  
  plot(muts_full$cor_sal, muts_full$af_cor_sal, ylim=c(-1,1), xlim=c(-1,1), xlab="SliM cor(p, Env2)", ylab="Subsample cor(p, Env2)", main=plotmain)
  abline(0,1)
  
  plot(muts_full$cor_temp, muts_full$cor_sal, main="SliM Correlation (AF, ENV) all samples", xlim=c(-1,1), ylim=c(-1,1))
  plot(muts_full$af_cor_temp, muts_full$af_cor_sal, main="Subsample Correlation (AF, ENV) all samples", xlim=c(-1,1), ylim=c(-1,1))
  
  
  
  muts_full$causal_temp <- "neutral-linked"
  muts_full$causal_temp[muts_full$causal & muts_full$mutTempEffect != 0] <- "causal"
  muts_full$causal_temp[muts_full$pos_pyslim > 500000] <- "neutral"
  
  ggplot(muts_full, aes(af_cor_temp, fill=causal_temp)) + 
    geom_density(alpha=0.5) + xlab("Cor(p, temp)") + xlim(-1,1) + ggtitle(plotmain) 
  
  
  muts_full$causal_sal <- "neutral-linked"
  muts_full$causal_sal[muts_full$causal & muts_full$mutSalEffect != 0 & info$Ntraits == 2] <- "causal"
  muts_full$causal_sal[muts_full$pos_pyslim > 500000] <- "neutral"
  
  ggplot(muts_full, aes(af_cor_sal, fill=causal_sal)) + 
    geom_density(alpha=0.5) + xlab("Core(p, Env2)") + xlim(-1,1) + ggtitle(plotmain)
  
  
  
  muts_full$pos_pyslim <- as.integer(muts_full$pos_pyslim)
  
  for (LG in 1:20){
    sites_low <- (LG-1)*50000+1
    sites_high <- (LG)*50000
    muts_full$LG[muts_full$pos_pyslim >= sites_low & muts_full$pos_pyslim <= sites_high] <- LG
  }
  
  muts_full$colors <- NA
  muts_full$colors[muts_full$LG %in% seq(1,9, by=2)] <- adjustcolor("goldenrod1", 0.1)
  muts_full$colors[muts_full$LG %in% seq(2,10, by=2)] <- adjustcolor("goldenrod3", 0.1)
  muts_full$colors[muts_full$LG %in% seq(11,19, by=2)] <- adjustcolor("grey50", 0.1)
  muts_full$colors[muts_full$LG %in% seq(12,20, by=2)] <- adjustcolor("grey70", 0.1)
  muts_full$colors <- as.character(muts_full$colors)
  
  muts_full$Va_temp <- muts_full$a_freq_subset*(1-muts_full$a_freq_subset)*muts_full$mutTempEffect^2
  
  muts_full$Va_temp_prop <- muts_full$Va_temp/sum(muts_full$Va_temp, na.rm=TRUE)
  
  hist(muts_full$Va_temp_prop[muts_full$causal_temp=="causal"], xlab="VA temp prop", main=plotmain, xlim=c(0,1), breaks=seq(0,1,0.01))
  muts_full$Va_temp_prop[is.na(muts_full$Va_temp_prop)] <- 0
  muts_full$Va_sal <- muts_full$a_freq_subset*(1-muts_full$a_freq_subset)*muts_full$mutSalEffect^2
  
  muts_full$Va_sal_prop <- muts_full$Va_sal/sum(muts_full$Va_sal, na.rm=TRUE)
  muts_full$Va_sal_prop[is.na(muts_full$Va_sal_prop)] <- 0
  hist(muts_full$Va_sal_prop[muts_full$causal_sal=="causal"], xlab= "VA Env2 prop.", main= plotmain, xlim=c(0,1), breaks=seq(0,1,0.01))
  
  plot(muts_full$Va_temp_prop[muts_full$causal==TRUE],
       muts_full$Va_sal_prop[muts_full$causal==TRUE],
       xlim=c(-0.01,1), ylim=c(-0.01,1), bty="l", pch=19, 
       xlab="Subset VA temp prop.",
       ylab="Subset VA Env2 prop.", main=plotmain)
  
  
  p1<- ggplot(aes(x=Va_temp_prop,y=abs(af_cor_temp), colour=causal_temp), data=muts_full) + 
    geom_point() +
    ggtitle(plotmain) + 
    theme_classic() + 
    scale_x_continuous(limits=c(0,1)) +
    scale_y_continuous(limits=c(0,1))
  ggExtra::ggMarginal(p1, type="violin",groupColour = TRUE, groupFill = TRUE)
  
  
  p<- ggplot(aes(x=Va_sal_prop,y=abs(af_cor_sal), colour=causal_sal), data=muts_full) + 
    geom_point() +
    ggtitle(plotmain) + 
    theme_classic() + 
    scale_x_continuous(limits=c(0,1)) +
    scale_y_continuous(limits=c(0,1))
  ggExtra::ggMarginal(p, type="violin",groupColour = TRUE, groupFill = TRUE)
  
  
  dev.off()
  
  
  ### Correlation stats ####
  Bonf_alpha <- (0.05/(nrow(muts_full)))
  
  muts_full$cor_temp_sig <- muts_full$af_cor_temp_P < Bonf_alpha
  muts_full$cor_sal_sig <- muts_full$af_cor_sal_P < Bonf_alpha
  
  num_causal_sig_temp_corr <- sum(muts_full$af_cor_temp_P[muts_full$causal_temp=="causal"]<Bonf_alpha)# number of causal loci that have significant Spearman's correlations with temperature after Bonferroni correction
  
  num_causal_sig_sal_corr<- sum(muts_full$af_cor_sal_P[muts_full$causal_sal=="causal"]<Bonf_alpha)# number of causal loci that have significant Spearman's correlations with salinity after Bonferroni correction
  
  num_notCausal_sig_temp_corr <- sum(muts_full$af_cor_temp_P[!muts_full$causal_temp=="causal"]<Bonf_alpha)# number of neutral loci that have significant  Spearman's correlations with temperature after Bonferroni correction
  num_notCausal_sig_sal_corr <- sum(muts_full$af_cor_sal_P[!muts_full$causal_sal=="causal"]<Bonf_alpha)# proportion of neutral loci that have significant Spearman's correlations with salinity after Bonferroni correction
  
  num_neut_sig_temp_corr <- sum(muts_full$af_cor_temp_P[!muts_full$causal_temp=="neutral"]<Bonf_alpha)# number of neutral loci that have significant  Spearman's correlations with temperature after Bonferroni correction
  num_neut_sig_sal_corr <- sum(muts_full$af_cor_sal_P[!muts_full$causal_sal=="neutral"]<Bonf_alpha)# proportion of neutral loci that have significant Spearman's correlations with salinity after Bonferroni correction
  
  cor_VA_temp_prop <- sum(muts_full$Va_temp_prop[which(muts_full$causal_temp=="causal" & (muts_full$af_cor_temp_P< Bonf_alpha) )]) 
  cor_VA_sal_prop <- sum(muts_full$Va_sal_prop[which(muts_full$causal_sal=="causal" & (muts_full$af_cor_sal_P< Bonf_alpha) )]) 
  
  cor_TPR_temp <-num_causal_sig_temp_corr/sum(muts_full$causal_temp=="causal")
  cor_FDR_allSNPs_temp <- num_notCausal_sig_temp_corr/sum(muts_full$af_cor_temp_P<Bonf_alpha) # all neutral / all outliers
  cor_FDR_neutSNPs_temp <- num_neut_sig_temp_corr/(num_neut_sig_temp_corr+num_causal_sig_temp_corr) # definitely neutral/ relevant outliers
  
  cor_TPR_sal <- num_causal_sig_sal_corr/sum(muts_full$causal_sal=="causal")
  cor_FDR_allSNPs_sal <- num_notCausal_sig_sal_corr/sum(muts_full$af_cor_sal_P<Bonf_alpha)
  cor_FDR_neutSNPs_sal <- num_neut_sig_sal_corr/(num_neut_sig_sal_corr+num_causal_sig_sal_corr)
  
  cor_AUCPR_temp_allSNPs<- (pr.curve(scores.class0 =-log10(muts_full$af_cor_temp_P[muts_full$causal_temp=="causal"]), scores.class1 = -log10(muts_full$af_cor_temp_P[!muts_full$causal_temp=="causal"])))$auc.integral
  # class 0 is the true positives, class 1 is the true negatives
  # it's important that the signal of class 0 is larger than class 1
  cor_AUCPR_temp_neutSNPs<- (pr.curve(scores.class0= -log10(muts_full$af_cor_temp_P[muts_full$causal_temp=="causal"]), scores.class1 = -log10(muts_full$af_cor_temp_P[muts_full$causal_temp=="neutral"])))$auc.integral
  
  cor_AUCPR_sal_allSNPs<- (pr.curve(scores.class0 =-log10(muts_full$af_cor_sal_P[muts_full$causal_sal=="causal"]), scores.class1 = -log10(muts_full$af_cor_sal_P[!muts_full$causal_sal=="causal"])))$auc.integral
  # class 0 is the true positives, class 1 is the true negatives
  # it's important that the signal of class 0 is larger than class 1
  cor_AUCPR_sal_neutSNPs<- (pr.curve(scores.class0= -log10(muts_full$af_cor_sal_P[muts_full$causal_sal=="causal"]), scores.class1 = -log10(muts_full$af_cor_sal_P[muts_full$causal_sal=="neutral"])))$auc.integral
  # PRROC: computing and visualizing precision-recall and receiver operating characteristic curves in R
  
  median_causal_temp_cor <- median(abs(muts_full$af_cor_temp[muts_full$causal_temp=="causal"]))#median abs(Spearman's correlation) between allele frequency and temperature for causal loci
  median_causal_sal_cor <- median(abs(muts_full$af_cor_sal[muts_full$causal_sal=="causal"])) #median Spearman's correlation  between allele frequency and salinity for causal loci
  
  median_neut_temp_cor <- median(abs(muts_full$af_cor_temp[!muts_full$causal_temp=="causal"])) # median Spearman's correlation  between allele frequency and temperature for neutral loci
  
  median_neut_sal_cor <-  median(abs(muts_full$af_cor_sal[!muts_full$causal_sal=="causal"]))# median Spearman's correlation  between allele frequency and salinity for neutral loci 
  
  
  ### Correlation Manhattan ####
  
  pdf(paste0(path,seed,"_pdf_8manhattan_cor.pdf"), width=10, height=5)
  
  cs <- "mako"
  begin_cs <- 0.9
  end_cs <- 0.2
  shape_causal <- 17
  outlier_color <- adjustcolor("brown1", 0.8)
  ymax=max(c(-log10(muts_full$af_cor_sal_P), -log10(muts_full$af_cor_temp_P), 10))
  
  a <- ggplot() + 
    geom_point(data=muts_full, aes(x=pos_pyslim, y = -log10(af_cor_temp_P)), color=muts_full$colors) +
    geom_point(data = muts_full[muts_full$causal_temp=="causal",], aes(x=pos_pyslim, y = -log10(af_cor_temp_P), 
                                                                       color=Va_temp_prop, size=Va_temp_prop), shape=shape_causal) + 
    scale_colour_viridis(option=cs, begin = begin_cs, end=end_cs, limits=c(0,1)) + 
    theme_classic() +
    geom_point(data = muts_full[muts_full$cor_temp_sig,], aes(x=pos_pyslim, y =-log10(af_cor_temp_P)),pch=23, col=outlier_color, size=3) + 
    ylim(0, ymax)  + ylab("-log10(P) Cor(af, temp)") +
    ggtitle(paste0(plotmain," temp")) + xlab("position") + labs(color="temp VA prop.", size="temp VA prop.")
  
  b <- ggplot() + 
    geom_point(data=muts_full, aes(x=pos_pyslim, y = -log10(af_cor_sal_P)), color=muts_full$colors) +
    geom_point(data = muts_full[muts_full$causal_sal=="causal",], aes(x=pos_pyslim, y = -log10(af_cor_sal_P), 
                                                                      color=Va_sal_prop, size=Va_sal_prop), shape=shape_causal) + 
    scale_colour_viridis(option=cs, begin = begin_cs, end=end_cs, limits=c(0,1)) + 
    theme_classic() +
    geom_point(data = muts_full[muts_full$cor_sal_sig,], aes(x=pos_pyslim, y =-log10(af_cor_sal_P)),pch=23, col=outlier_color, size=3) + 
    ylim(0, ymax)  + ylab("-log10(P) Cor(af, Env2)") +
    ggtitle(paste0(plotmain," temp")) + xlab("position") + labs(color="Env2 VA prop.", size="Env2 VA prop.")
  
  grid.arrange(a, b, nrow=2)
  
  dev.off()
  
### PCA ####
  lfmmfile <- paste0(path, seed, "_genotypes.lfmm")
  write.lfmm(t(G_full_subset), lfmmfile)
  
  pc = pca(lfmmfile, 30, scale = TRUE)
  
  subset_indPhen_df$PC1 <- pc$projections[,1]
  subset_indPhen_df$PC2 <- pc$projections[,2]
  subset_indPhen_df$PC3 <- pc$projections[,3]
  tw = tracy.widom(pc)
  plot(tw$percentage)
  a <- tw$percentage[1:20]
  b <- tw$percentage[2:21]
  K <- max(which(a > b*1.5))
  K
  
  pdf(paste0(path,seed,"_pdf_3pca.pdf"), width=5, height=5)
  
  plot(pc$sdev[1:10]/sum(pc$sdev), bty="l", ylab="Prop Var of PC axis", main=paste0(plotmain, "; K=", K))
  
  ggplot(subset_indPhen_df) + geom_point(aes(x=PC1, y=PC2,size=sal_opt)) + geom_point(aes(x=PC1, y=PC2, color=temp_opt)) + theme_classic() + scale_colour_gradient2(high=rgb(1,0.4,0.2), low="cornflowerblue", mid=rgb(0.8,0.8,0.7), name="Temp") + labs(size="Env2") + ggtitle(paste0(plotmain, "; K=", K))
  
  ggplot(subset_indPhen_df) + geom_point(aes(x=PC1, y=PC3,size=sal_opt)) + geom_point(aes(x=PC1, y=PC3, color=temp_opt)) + theme_classic() + scale_colour_gradient2(high=rgb(1,0.4,0.2), low="cornflowerblue", mid=rgb(0.8,0.8,0.7), name="Temp") + labs(size="Env2") + ggtitle(paste0(plotmain,  "; K=", K))
  
  
  dev.off()
  
# LFMM ####
  
  mod_temp <- lfmm2(input = t(G_full_subset), env = subset_indPhen_df$temp_opt, K = K)
  pv_temp <- lfmm2.test(object = mod_temp, input = t(G_full_subset),
                        env = subset_indPhen_df$temp_opt,
                        linear = TRUE)
  
  muts_full$LEA3.2_lfmm2_mlog10P_tempenv <- -log10(as.numeric(pv_temp$pvalues))
  muts_full$LEA3.2_lfmm2_mlog10P_tempenv_sig <- qvalue(as.numeric(pv_temp$pvalues))$qvalue < 0.05
  mod_sal <- lfmm2(input = t(G_full_subset), env = subset_indPhen_df$sal_opt, K = K)
  
  pv_sal <- lfmm2.test(object = mod_sal, input = t(G_full_subset),
                       env = subset_indPhen_df$sal_opt,
                       linear = TRUE)
  
  
  muts_full$LEA3.2_lfmm2_mlog10P_salenv <- -log10(as.numeric(pv_sal$pvalues))
  muts_full$LEA3.2_lfmm2_mlog10P_salenv_sig <- qvalue(as.numeric(pv_sal$pvalues))$qvalue < 0.05
  
  LEA3.2_lfmm2_Va_temp_prop <- sum(muts_full$Va_temp[which(muts_full$LEA3.2_lfmm2_mlog10P_tempenv_sig)], na.rm=TRUE)/sum(muts_full$Va_temp, na.rm=TRUE) 
  LEA3.2_lfmm2_Va_sal_prop <- sum(muts_full$Va_sal[which(muts_full$LEA3.2_lfmm2_mlog10P_salenv_sig)], na.rm=TRUE)/sum(muts_full$Va_sal, na.rm=TRUE) 
  
  # True positive rate -proportion of causal loci that are outliers
  LEA3.2_lfmm2_TPR_temp <- sum(muts_full$LEA3.2_lfmm2_mlog10P_tempenv_sig & muts_full$causal_temp=="causal")/sum(muts_full$causal_temp=="causal")
  
  LEA3.2_lfmm2_TPR_sal <- sum(muts_full$LEA3.2_lfmm2_mlog10P_salenv_sig & muts_full$causal_sal=="causal")/sum(muts_full$causal_sal=="causal")
  
  # False discovery rate- proportion of outliers that are neutral
  LEA3.2_lfmm2_FDR_allSNPs_temp <-(sum(!muts_full$causal_temp=="causal" & muts_full$LEA3.2_lfmm2_mlog10P_tempenv_sig, na.rm=TRUE)/sum(muts_full$LEA3.2_lfmm2_mlog10P_tempenv_sig, na.rm=TRUE))
  LEA3.2_lfmm2_FDR_allSNPs_sal <- (sum(!muts_full$causal_sal=="causal" & muts_full$LEA3.2_lfmm2_mlog10P_salenv_sig, na.rm=TRUE)/sum(muts_full$LEA3.2_lfmm2_mlog10P_salenv_sig, na.rm=TRUE))
    # non-causal (all neutral) / all outliers
  
  LEA3.2_lfmm2_num_causal_sig_temp <- sum(muts_full$causal_temp=="causal" & muts_full$LEA3.2_lfmm2_mlog10P_tempenv_sig, na.rm=TRUE)
  LEA3.2_lfmm2_num_neut_sig_temp <-  sum(muts_full$causal_temp=="neutral" & muts_full$LEA3.2_lfmm2_mlog10P_tempenv_sig, na.rm=TRUE)
  
  LEA3.2_lfmm2_num_causal_sig_sal <- sum(muts_full$causal_sal=="causal" & muts_full$LEA3.2_lfmm2_mlog10P_salenv_sig, na.rm=TRUE)
  LEA3.2_lfmm2_num_neut_sig_sal <-  sum(muts_full$causal_sal=="neutral" & muts_full$LEA3.2_lfmm2_mlog10P_salenv_sig, na.rm=TRUE)
  
  LEA3.2_lfmm2_FDR_neutSNPs_temp <-LEA3.2_lfmm2_num_neut_sig_temp/(LEA3.2_lfmm2_num_neut_sig_temp+LEA3.2_lfmm2_num_causal_sig_temp )
  LEA3.2_lfmm2_FDR_neutSNPs_sal <- LEA3.2_lfmm2_num_neut_sig_sal/(LEA3.2_lfmm2_num_neut_sig_sal+LEA3.2_lfmm2_num_causal_sig_sal )
    # neutral / relevant outliers
  
  LEA3.2_lfmm2_AUCPR_temp_allSNPs<- (pr.curve(scores.class0 =muts_full$LEA3.2_lfmm2_mlog10P_salenv[muts_full$causal_temp=="causal"], scores.class1 = muts_full$LEA3.2_lfmm2_mlog10P_salenv[!muts_full$causal_temp=="causal"]))$auc.integral
  LEA3.2_lfmm2_AUCPR_temp_neutSNPs<- (pr.curve(scores.class0 =muts_full$LEA3.2_lfmm2_mlog10P_salenv[muts_full$causal_temp=="causal"], scores.class1 = muts_full$LEA3.2_lfmm2_mlog10P_salenv[muts_full$causal_temp=="neutral"]))$auc.integral
  LEA3.2_lfmm2_AUCPR_sal_allSNPs<- (pr.curve(scores.class0 =muts_full$LEA3.2_lfmm2_mlog10P_salenv[muts_full$causal_sal=="causal"], scores.class1 = muts_full$LEA3.2_lfmm2_mlog10P_salenv[!muts_full$causal_sal=="causal"]))$auc.integral
  LEA3.2_lfmm2_AUCPR_sal_neutSNPs<- (pr.curve(scores.class0 =muts_full$LEA3.2_lfmm2_mlog10P_salenv[muts_full$causal_sal=="causal"], scores.class1 = muts_full$LEA3.2_lfmm2_mlog10P_salenv[muts_full$causal_sal=="neutral"]))$auc.integral
  # PRROC: computing and visualizing precision-recall and receiver operating characteristic curves in R
  # class 0 is the true positives, class 1 is the true negatives
  # it's important that the signal of class 0 is larger than class 1
  
  
## Correlate causal loci AF with PCA and LFMM pop structure ####
  muts_full$structure_cor_G_PC1 <- muts_full$structure_cor_G_LFMM_U1_modtemp<- muts_full$structure_cor_G_LFMM_U1_modsal<- NA
  for (i in 1:nrow(muts_full)){
    muts_full$structure_cor_G_PC1[i] <- cor(G_full_subset[i,], subset_indPhen_df$PC1, method = "kendall")
    muts_full$structure_cor_G_LFMM_U1_modtemp[i] <- cor(G_full_subset[i,],mod_temp@U[,1], method = "kendall")
    muts_full$structure_cor_G_LFMM_U1_modsal[i] <- cor(G_full_subset[i,],mod_sal@U[,1], method = "kendall")  
  }

  #VIZ:
  #sig by different methods vs. corr locus with structure
  #histogram of causal vs. neutral loci correlation with structure
  
  # cor plot of the different structure correlations
    par(mfrow=c(2,1))
    subset_indPhen_df$LFMM_U1_temp <- mod_temp@U[,1]
    subset_indPhen_df$LFMM_U1_sal <- mod_sal@U[,1]
    #plot(subset_indPhen_df$PC1, mod_temp@U[,1])
    #plot(subset_indPhen_df$PC1, mod_sal@U[,1])
    
    # PC vs. Env
    k <- ggplot(subset_indPhen_df) + geom_point(aes(x=PC1, y=temp_opt, size=sal_opt)) + geom_point(aes(x=PC1, y=temp_opt, color=temp_opt)) + ylab("Deme temperature") + theme_classic() + scale_colour_gradient2(high=rgb(1,0.4,0.2), low="cornflowerblue", mid=rgb(0.8,0.8,0.7), name="Temp") + labs(size="Env2") + ggtitle(paste0(plotmain))
    l <- ggplot(subset_indPhen_df) + geom_point(aes(x=PC1, y=sal_opt, size=sal_opt)) + geom_point(aes(x=PC1, y=sal_opt, color=temp_opt)) + ylab("Deme Env2") + theme_classic() + scale_colour_gradient2(high=rgb(1,0.4,0.2), low="cornflowerblue", mid=rgb(0.8,0.8,0.7), name="Temp") + labs(size="Env2") + ggtitle(paste0(plotmain))
    
    # PC vs LFMM latent factors
    m <- ggplot(subset_indPhen_df) + geom_point(aes(x=PC1, y=LFMM_U1_temp,size=sal_opt)) + geom_point(aes(x=PC1, y=LFMM_U1_temp, color=temp_opt)) + ylab("Latent factor 1 LFMM temp model") + theme_classic() + scale_colour_gradient2(high=rgb(1,0.4,0.2), low="cornflowerblue", mid=rgb(0.8,0.8,0.7), name="Temp") + labs(size="Env2") + ggtitle(paste0(plotmain))
    n <- ggplot(subset_indPhen_df) + geom_point(aes(x=PC1, y=LFMM_U1_sal,size=sal_opt)) + geom_point(aes(x=PC1, y=LFMM_U1_sal, color=temp_opt)) + ylab("Latent factor 1 LFMM Env2 model") + theme_classic() + scale_colour_gradient2(high=rgb(1,0.4,0.2), low="cornflowerblue", mid=rgb(0.8,0.8,0.7), name="Temp") + labs(size="Env2") + ggtitle(paste0(plotmain))
    o <- ggplot(subset_indPhen_df) + geom_point(aes(x=PC2, y=LFMM_U1_temp,size=sal_opt)) + geom_point(aes(x=PC2, y=LFMM_U1_temp, color=temp_opt)) + ylab("Latent factor 1 LFMM temp model") + theme_classic() + scale_colour_gradient2(high=rgb(1,0.4,0.2), low="cornflowerblue", mid=rgb(0.8,0.8,0.7), name="Temp") + labs(size="Env2") + ggtitle(paste0(plotmain))
    p <- ggplot(subset_indPhen_df) + geom_point(aes(x=PC2, y=LFMM_U1_sal,size=sal_opt)) + geom_point(aes(x=PC2, y=LFMM_U1_sal, color=temp_opt)) + ylab("Latent factor 1 LFMM Env2 model") + theme_classic() + scale_colour_gradient2(high=rgb(1,0.4,0.2), low="cornflowerblue", mid=rgb(0.8,0.8,0.7), name="Temp") + labs(size="Env2") + ggtitle(paste0(plotmain))
    grid.arrange(k, l, m, n, o, p, nrow=3)
  
  # Plot Deme temperature of individual vs. PC1 loading of individual
  # Plot cor(Genotype, structure PC1) vs. -log10 P cor(af, temp)
  # Plot cor(Genotype, structure PC1) vs. -log10 P lfmm temp model

    # PC1 plots
    # Temp corrected
    r <- ggplot() + 
      geom_point(data=muts_full, aes(x=abs(structure_cor_G_PC1), y = LEA3.2_lfmm2_mlog10P_tempenv), color=muts_full$colors) +
      geom_point(data = muts_full[muts_full$causal_temp=="causal",], aes(x=abs(structure_cor_G_PC1), y = LEA3.2_lfmm2_mlog10P_tempenv, 
                                                                         color=Va_temp_prop, size=Va_temp_prop), shape=shape_causal) + 
      scale_colour_viridis(option=cs, begin = begin_cs, end=end_cs, limits=c(0,1)) + 
      theme_classic() +
      geom_point(data = muts_full[muts_full$LEA3.2_lfmm2_mlog10P_tempenv_sig,], aes(x=abs(structure_cor_G_PC1), y =LEA3.2_lfmm2_mlog10P_tempenv),pch=23, col=outlier_color, size=3) + 
      ylim(0, ymax)  + ylab("-log10(P) LFMM (Genotype, temp)") +
      ggtitle("Temperature - Structure corrected") + xlab("Abs(Cor(Genotype, Structure PC1))") + labs(color="temp VA prop.", size="temp VA prop.")
    
    # Temp uncorreected
    s <- ggplot() + 
      geom_point(data=muts_full, aes(x=abs(structure_cor_G_PC1), y = -log10(af_cor_temp_P)), color=muts_full$colors) + 
      geom_point(data = muts_full[muts_full$causal_temp=="causal",], aes(x=abs(structure_cor_G_PC1), y = -log10(af_cor_temp_P), 
                                                                         color=Va_temp_prop, size=Va_temp_prop), shape=shape_causal) + 
      scale_colour_viridis(option=cs, begin = begin_cs, end=end_cs, limits=c(0,1)) + 
      theme_classic()  +
      geom_point(data = muts_full[muts_full$cor_temp_sig,], aes(x=abs(structure_cor_G_PC1), y = -log10(af_cor_temp_P)),pch=23, col=adjustcolor("darkorange",0.5), size=3) + 
      ylim(0, ymax)  + ylab("-log10(P) cor(Genotype, temp)") +
      ggtitle("Temperature - Structure not corrected") + xlab("Abs(Cor(Genotype, Structure PC1))") + labs(color="temp VA prop.", size="temp VA prop.")
   
    # Sal corrected
    rsal <- ggplot() + 
      geom_point(data=muts_full, aes(x=abs(structure_cor_G_PC1), y = LEA3.2_lfmm2_mlog10P_salenv), color=muts_full$colors) +
      geom_point(data = muts_full[muts_full$causal_sal=="causal",], aes(x=abs(structure_cor_G_PC1), y = LEA3.2_lfmm2_mlog10P_salenv, 
                                                                         color=Va_sal_prop, size=Va_sal_prop), shape=shape_causal) + 
      scale_colour_viridis(option=cs, begin = begin_cs, end=end_cs, limits=c(0,1)) + 
      theme_classic() +
      geom_point(data = muts_full[muts_full$LEA3.2_lfmm2_mlog10P_salenv_sig,], aes(x=abs(structure_cor_G_PC1), y =LEA3.2_lfmm2_mlog10P_salenv),pch=23, col=outlier_color, size=3) + 
      ylim(0, ymax)  + ylab("-log10(P) LFMM (Genotype, Env2)") +
      ggtitle("Env2 - Structure corrected") + xlab("Abs(Cor(Genotype, Structure PC1))") + labs(color="Env2 VA prop.", size="Env2 VA prop.")
    
    # Sal uncorrected
    ssal <- ggplot() + 
      geom_point(data=muts_full, aes(x=abs(structure_cor_G_PC1), y = -log10(af_cor_sal_P)), color=muts_full$colors) + 
      geom_point(data = muts_full[muts_full$causal_sal=="causal",], aes(x=abs(structure_cor_G_PC1), y = -log10(af_cor_sal_P), 
                                                                         color=Va_sal_prop, size=Va_sal_prop), shape=shape_causal) + 
      scale_colour_viridis(option=cs, begin = begin_cs, end=end_cs, limits=c(0,1)) + 
      theme_classic()  +
       geom_point(data = muts_full[muts_full$cor_sal_sig,], aes(x=abs(structure_cor_G_PC1), y =-log10(af_cor_sal_P)),pch=23, col=adjustcolor("darkorange",0.5), size=3) + 
      ylim(0, ymax)  + ylab("-log10(P) cor(Genotype, Env2)") +
      ggtitle("Env2 - Structure not corrected") + xlab("Abs(Cor(Genotype, Structure PC1))") + labs(color="Env2 VA prop.", size="Env2 VA prop.")
    
    
    grid.arrange(k, l, s, ssal, r, rsal, nrow=3)
  
    #cor.test(abs(muts_full$structure_cor_G_PC1), muts_full$LEA3.2_lfmm2_mlog10P_tempenv, method="pearson")
    
    # q <- ggplot() + 
    #   geom_point(data=muts_full, aes(x=abs(structure_cor_G_LFMM_U1_modtemp), y = LEA3.2_lfmm2_mlog10P_tempenv), color=muts_full$colors) +
    #   geom_point(data = muts_full[muts_full$causal_temp=="causal",], aes(x=abs(structure_cor_G_LFMM_U1_modtemp), y = LEA3.2_lfmm2_mlog10P_tempenv, 
    #                                                                      color=Va_temp_prop, size=Va_temp_prop), shape=shape_causal) + 
    #   scale_colour_viridis(option=cs, begin = begin_cs, end=end_cs, limits=c(0,1)) + 
    #   theme_classic() +
    #   geom_point(data = muts_full[muts_full$LEA3.2_lfmm2_mlog10P_tempenv_sig,], aes(x=abs(structure_cor_G_LFMM_U1_modtemp), y =LEA3.2_lfmm2_mlog10P_tempenv),pch=23, col=outlier_color, size=3) + 
    #   ylim(0, ymax)  + ylab("-log10(P) LFMM temp") +
    #   ggtitle(paste0(plotmain," temp")) + xlab("Abs(Cor(Genotype, Structure LFMM U1 temp model))") + labs(color="temp VA prop.", size="temp VA prop.")
    # 
    # qaf <- ggplot() + 
    #   geom_point(data=muts_full, aes(x=abs(structure_cor_G_LFMM_U1_modtemp), y = -log10(af_cor_temp_P)), color=muts_full$colors) +
    #   geom_point(data = muts_full[muts_full$causal_temp=="causal",], aes(x=abs(structure_cor_G_LFMM_U1_modtemp), y = -log10(af_cor_temp_P), 
    #                                                                      color=Va_temp_prop, size=Va_temp_prop), shape=shape_causal) + 
    #   scale_colour_viridis(option=cs, begin = begin_cs, end=end_cs, limits=c(0,1)) + 
    #   theme_classic() +
    #   #geom_point(data = muts_full[muts_full$LEA3.2_lfmm2_mlog10P_tempenv_sig,], aes(x=abs(structure_cor_G_LFMM_U1_modtemp), y =-log10(af_cor_temp_P),pch=23, col=outlier_color, size=3) + 
    #   ylim(0, ymax)  + ylab("-log10(P) LFMM temp") +
    #   ggtitle(paste0(plotmain," temp")) + xlab("Abs(Cor(Genotype, Structure LFMM U1 temp model))") + labs(color="temp VA prop.", size="temp VA prop.")
    # 
    # grid.arrange(qaf, q, nrow=2)
    # T ODO FIX Y AXES
    
    # t <- ggplot() + 
    #   geom_point(data=muts_full, aes(y=LEA3.2_lfmm2_mlog10P_tempenv, x = -log10(af_cor_temp_P)), color=muts_full$colors) +
    #   geom_point(data = muts_full[muts_full$causal_temp=="causal",], aes(y=LEA3.2_lfmm2_mlog10P_tempenv, x = -log10(af_cor_temp_P), 
    #                                                                      color=Va_temp_prop, size=Va_temp_prop), shape=shape_causal) + 
    #   scale_colour_viridis(option=cs, begin = begin_cs, end=end_cs, limits=c(0,1)) + 
    #   theme_classic() +
    #   ylim(0, ymax)  + ylab("-log10(P) LFMM temp") +
    #   ggtitle(paste0(plotmain," temp")) + xlab("-log10 P cor(af, temp)") + labs(color="temp VA prop.", size="temp VA prop.")
    # 
  
  ### LFMM Manhattan ####
  
  pdf(paste0(path,seed,"_pdf_4manhattan_LFMM.pdf"), width=10, height=5)
  
  cs <- "mako"
  begin_cs <- 0.9
  end_cs <- 0.2
  shape_causal <- 17
  outlier_color <- adjustcolor("brown1", 0.8)
  ymax=max(c(muts_full$LEA3.2_lfmm2_mlog10P_tempenv, muts_full$LEA3.2_lfmm2_mlog10P_salenv, 10))
  
  a <- ggplot() + 
    geom_point(data=muts_full, aes(x=pos_pyslim, y = LEA3.2_lfmm2_mlog10P_tempenv), color=muts_full$colors) +
    geom_point(data = muts_full[muts_full$causal_temp=="causal",], aes(x=pos_pyslim, y = LEA3.2_lfmm2_mlog10P_tempenv, 
                                                                       color=Va_temp_prop, size=Va_temp_prop), shape=shape_causal) + 
    scale_colour_viridis(option=cs, begin = begin_cs, end=end_cs, limits=c(0,1)) + 
    theme_classic() +
    geom_point(data = muts_full[muts_full$LEA3.2_lfmm2_mlog10P_tempenv_sig,], aes(x=pos_pyslim, y =LEA3.2_lfmm2_mlog10P_tempenv),pch=23, col=outlier_color, size=3) + 
    ylim(0, ymax)  + ylab("-log10(P) LFMM temp") +
    ggtitle(paste0(plotmain," temp")) + xlab("position") + labs(color="temp VA prop.", size="temp VA prop.")
  
  b <- ggplot() + 
    geom_point(data=muts_full, aes(x=pos_pyslim, y = LEA3.2_lfmm2_mlog10P_salenv), color=muts_full$colors) +
    geom_point(data = muts_full[muts_full$causal_sal=="causal",], aes(x=pos_pyslim, y = LEA3.2_lfmm2_mlog10P_salenv, 
                                                                      color=Va_sal_prop, size=Va_sal_prop), shape=shape_causal) + 
    scale_colour_viridis(option=cs, begin = begin_cs, end=end_cs, limits=c(0,1)) + 
    theme_classic() +
    geom_point(data = muts_full[muts_full$LEA3.2_lfmm2_mlog10P_salenv_sig,], aes(x=pos_pyslim, y =LEA3.2_lfmm2_mlog10P_salenv),pch=23, col=outlier_color, size=3) + 
    ylim(0, ymax)  + ylab("-log10(P) LFMM salinity") + ggtitle(paste0(plotmain," Env2")) + xlab("position") + labs(color="Env2 VA prop.", size="Env2 VA prop.")
  
  grid.arrange(a, b, nrow=2)
  
  dev.off()
  
### RDA vegan ####
  
  sal <- subset_indPhen_df$sal_opt
  temp <- subset_indPhen_df$temp_opt
  rdaout <- rda(t(G_full_subset)~ sal + temp)
  
  str(rdaout)
  scores <- scores(rdaout, choices=1:4)
  loci.sc <- scores$species
  ind.sc <- scores$sites
  
  subset_indPhen_df$RDA1 <- ind.sc[,1]
  subset_indPhen_df$RDA2 <- ind.sc[,2]
  
  muts_full$RDA1_score <- loci.sc[,1]
  muts_full$RDA2_score <- loci.sc[,2]
  
  rdadapt<-function(rda,K) {
    zscores<-rda$CCA$v[,1:as.numeric(K)]
    resscale <- apply(zscores, 2, scale)
    resmaha <- covRob(resscale, distance = TRUE, na.action= na.omit, estim="pairwiseGK")$dis 
    lambda <- median(resmaha)/qchisq(0.5,df=K)
    reschi2test <- pchisq(resmaha/lambda,K,lower.tail=FALSE)
    qval <- qvalue(reschi2test)
    q.values_rdadapt<-qval$qvalues
    return(data.frame(p.values=reschi2test, q.values=q.values_rdadapt))
  }
  
# RDA outliers and error rates####
  
  ps <- rdadapt(rdaout, 2)
  
  muts_full$RDA_mlog10P <- -log10(ps$p.values)
  muts_full$RDA_mlog10P_sig <- ps$q.values<0.001
  
  
  
  # Proportion VA
  (RDA_Va_temp_prop <- sum(muts_full$Va_temp_prop[muts_full$RDA_mlog10P_sig], na.rm=TRUE))
  
  (RDA_Va_sal_prop <- sum(muts_full$Va_sal_prop[muts_full$RDA_mlog10P_sig], na.rm=TRUE))
  
  (RDA_TPR <- sum(muts_full$RDA_mlog10P_sig & muts_full$causal)/sum(muts_full$causal))
  (RDA_FDR_allSNPs <- sum(muts_full$RDA_mlog10P_sig & !muts_full$causal)/sum(muts_full$RDA_mlog10P_sig)) #non-causal outliers / (all outliers)
  
  num_RDA_sig_causal <- sum(muts_full$RDA_mlog10P_sig & muts_full$causal)
  num_RDA_sig_neutral <- sum(muts_full$RDA_mlog10P_sig & muts_full$causal_temp=="neutral")
  RDA_FDR_neutSNPs <- num_RDA_sig_neutral/(num_RDA_sig_neutral + num_RDA_sig_causal) #neutral outliers / (all outliers)
  
  RDA_AUCPR_allSNPs<- (pr.curve(scores.class0 = muts_full$RDA_mlog10P[muts_full$causal], scores.class1 = muts_full$RDA_mlog10P[!muts_full$causal]))$auc.integral
  # class 0 is the true positives, class 1 is the true negatives
  # it's important that the signal of class 0 is larger than class 1
  RDA_AUCPR_neutSNPs<- (pr.curve(scores.class0= muts_full$RDA_mlog10P[muts_full$causal], scores.class1 = muts_full$RDA_mlog10P[muts_full$causal_temp=="neutral"]))$auc.integral
    # only neutral SNPs on 2nd half of genome not affected by selection
  
# RDA outlier and individual plots ####
  muts_scale <- c(max(abs(muts_full$RDA1_score))*1.3, max(abs(muts_full$RDA2_score))*1.3)
  # for nicer plot
  sal_arrow_muts = rdaout$CCA$biplot[which(rownames(rdaout$CCA$biplot)=="sal"),]*muts_scale
  temp_arrow_muts = rdaout$CCA$biplot[which(rownames(rdaout$CCA$biplot)=="temp"),]*muts_scale
  
  arrows <- data.frame(x=rep(0,2), y = rep(0, 2),
                       dx = c(sal_arrow_muts[1], temp_arrow_muts[1]),
                       dy = c(sal_arrow_muts[2], temp_arrow_muts[2]), name = c("Env2", "Temp"))
  
  pdf(paste0(path,seed,"_pdf_5RDA.pdf"), width=5, height=5)
  
  a<- screeplot(rdaout)
  str(a)
  a$y # save this it's the eigenvalues
  
  ## Mutation RDA
  
  # ggplot() + 
  #   geom_point(data=muts_full, aes(x=RDA1_score, y = RDA2_score), color=adjustcolor("grey80", 0.2)) +
  #   geom_point(data = muts_full[muts_full$causal_temp=="causal",], 
  #              aes(x=RDA1_score, y = RDA2_score, color=Va_temp_prop, size=Va_temp_prop), shape=shape_causal) +
  #   scale_colour_viridis(option=cs, begin = begin_cs, end=end_cs, limits=c(0,1)) +  
  #   theme_classic() +  geom_point(data = muts_full[muts_full$RDA_mlog10P_sig,], aes(x=RDA1_score, y = RDA2_score),
  #                                 pch=23, col=outlier_color, size=3) +  xlab("RDA 1") + ylab("RDA 2") + 
  #   ggtitle(paste0(plotmain, " mutations - temp")) + 
  #   geom_segment(data=arrows, aes(x=x, y=y, xend=dx, yend=dy),arrow=arrow(length = unit(0.2,"cm"))) + geom_text(data=arrows, aes(x=dx, y=dy, label=name), hjust="right", vjust="bottom")+ labs(color="temp VA prop.", size="temp VA prop.")
  # 
  
  lims <- max(abs(c(muts_full$mutSalEffect,muts_full$mutTempEffect)), na.rm=TRUE)
  
  ggplot() + 
    geom_point(data=muts_full, aes(x=RDA1_score, y = RDA2_score), color=adjustcolor("grey80", 0.2)) +
    geom_point(data = muts_full[muts_full$causal_temp=="causal",], 
               aes(x=RDA1_score, y = RDA2_score, color=mutTempEffect, size=Va_temp_prop), shape=shape_causal) +
    scale_colour_viridis(option="turbo", #begin = begin_cs, end=end_cs, 
                         limits=c(-lims, lims)
                         )+  
    theme_classic() +  geom_point(data = muts_full[muts_full$RDA_mlog10P_sig,], aes(x=RDA1_score, y = RDA2_score),
                                  pch=23, col=adjustcolor("black",0.5), size=3) +  xlab("RDA 1") + ylab("RDA 2") + 
    ggtitle(paste0(plotmain, " mutations - temp")) + 
    geom_segment(data=arrows, aes(x=x, y=y, xend=dx, yend=dy),arrow=arrow(length = unit(0.2,"cm"))) + 
    geom_text(data=arrows, aes(x=dx, y=dy, label=name), hjust="right", vjust="bottom")+ labs(color="mut Temp effect", size="temp VA prop.")
  
  
  ggplot() + 
    geom_point(data=muts_full, aes(x=RDA1_score, y = RDA2_score), color=adjustcolor("grey80", 0.2)) +
    geom_point(data = muts_full[muts_full$causal_sal=="causal",], 
               aes(x=RDA1_score, y = RDA2_score, color=mutSalEffect, size=Va_sal_prop), shape=shape_causal) +
    scale_colour_viridis(option="turbo", #begin = begin_cs, end=end_cs, 
                         limits=c(-lims, lims)
    )+  
    theme_classic() +  geom_point(data = muts_full[muts_full$RDA_mlog10P_sig,], aes(x=RDA1_score, y = RDA2_score),
                                  pch=23, col=adjustcolor("black",0.5), size=3) +  xlab("RDA 1") + ylab("RDA 2") + 
    ggtitle(paste0(plotmain, " mutations - Env2")) + 
    geom_segment(data=arrows, aes(x=x, y=y, xend=dx, yend=dy),arrow=arrow(length = unit(0.2,"cm"))) + 
    geom_text(data=arrows, aes(x=dx, y=dy, label=name), hjust="right", vjust="bottom")+ labs(color="mut Env2 effect", size="Env2 VA prop.")
  
  
## RDA individual plot
  ind_scale <- c(max(abs(subset_indPhen_df$RDA1))*1.3, max(abs(subset_indPhen_df$RDA2))*1.3)
  # for nicer plot
  sal_arrow_muts = rdaout$CCA$biplot[which(rownames(rdaout$CCA$biplot)=="sal"),]*ind_scale
  temp_arrow_muts = rdaout$CCA$biplot[which(rownames(rdaout$CCA$biplot)=="temp"),]*ind_scale
  
  arrows_ind <- data.frame(x=rep(0,2), y = rep(0, 2),
                           dx = c(sal_arrow_muts[1], temp_arrow_muts[1]),
                           dy = c(sal_arrow_muts[2], temp_arrow_muts[2]), name = c("Env2", "Temp"))
  
  ggplot(subset_indPhen_df) + 
    geom_point(aes(x=RDA1, y = RDA2, size=sal_opt)) + geom_point(aes(x=RDA1, y=RDA2, color=temp_opt)) + theme_classic() + 
    scale_colour_gradient2(high=rgb(1,0.4,0.2), low="cornflowerblue", mid=rgb(0.8,0.8,0.7), name="Temp") + labs(size="Env2") + 
    ggtitle(paste0(plotmain,"; Individs.; N traits = ", info$Ntraits)) + 
    geom_segment(data=arrows_ind, aes(x=x, y=y, xend=dx, yend=dy),arrow=arrow(length = unit(0.2,"cm"))) + 
    geom_text(data=arrows_ind, aes(x=dx, y=dy, label=name), hjust="right", vjust="bottom")
  
  dev.off()
  
#### RDA manahattan ####
  ymax = max(c(muts_full$RDA_mlog10P, 10))
  pdf(paste0(path,seed,"_pdf_6manhattan_RDA.pdf"), width=10, height=5)
  
  a<- ggplot() + 
    geom_point(data=muts_full, aes(x=pos_pyslim, y = RDA_mlog10P), color=muts_full$colors) +
    geom_point(data = muts_full[muts_full$causal_temp=="causal",], 
               aes(x=pos_pyslim, y = RDA_mlog10P, color=Va_temp_prop, size=Va_temp_prop), shape=shape_causal) + 
    scale_colour_viridis(option=cs, begin = begin_cs, end=end_cs, limits=c(0,1)) + 
    theme_classic() +
    geom_point(data = muts_full[muts_full$RDA_mlog10P_sig,], aes(x=pos_pyslim, y =RDA_mlog10P),pch=23, col=outlier_color, size=3) + 
    ylim(0, ymax)  + ylab("-log10(P) RDA") + ggtitle(paste0(plotmain," temp")) + xlab("position") + labs(color="Temp VA prop.", size="Temp VA prop.")
  
  b<- ggplot() + 
    geom_point(data=muts_full, aes(x=pos_pyslim, y = RDA_mlog10P), color=muts_full$colors) +
    geom_point(data = muts_full[muts_full$causal_sal=="causal",], aes(x=pos_pyslim, y = RDA_mlog10P, color=Va_sal_prop, size=Va_sal_prop), shape=shape_causal) + 
    scale_colour_viridis(option=cs, begin =begin_cs, end=end_cs, limits=c(0,1)) + 
    theme_classic() +
    geom_point(data = muts_full[muts_full$RDA_mlog10P_sig,], aes(x=pos_pyslim, y =RDA_mlog10P),pch=23, col=outlier_color, size=3) + 
    ylim(0, ymax)  + ylab("-log10(P) RDA") + ggtitle(paste0(plotmain," Env2")) + xlab("position") + labs(color="Env2 VA prop.", size="Env2 VA prop.")
  
  grid.arrange(a, b, nrow=2)
  dev.off()
  
#### write RDA table ####
  write.table(data.frame(seed=as.character(seed), summary(rdaout)$biplot),paste0(path,seed,"_RDAloadings.txt")) 
  
### RDA predict phenotype from random loci ####
  nmax <- min(c(20000, num_causal+num_neut))
  nloci <- c(10, 50, 100, 500, nmax)
  Random_cor_salpredict_salphen <- rep(NA, length(nloci))
  Random_cor_temppredict_tempphen <- rep(NA, length(nloci))
  
  for (i in 1:length(nloci)){
    G_random_out <- G_full_subset[sample(1:nrow(muts_full), nloci[i], replace=FALSE),]
    rdaout_rand <- rda(t(G_random_out)~ subset_indPhen_df$sal_opt + subset_indPhen_df$temp_opt)
    scores <- scores(rdaout_rand, choices=1:4)
    loci.sc <- scores$species
    ind.sc <- scores$sites
    temp_pred <- ind.sc[,1]*eigenvals(rdaout_rand)[1]*summary(rdaout_rand)$biplot[2,1] + ind.sc[,2]*eigenvals(rdaout_rand)[2]*summary(rdaout_rand)$biplot[2,2]
    sal_pred <- ind.sc[,1]*eigenvals(rdaout_rand)[1]*summary(rdaout_rand)$biplot[1,1] + ind.sc[,2]*eigenvals(rdaout_rand)[2]*summary(rdaout_rand)$biplot[1,2]
    
    Random_cor_temppredict_tempphen[i] <- cor(scale(subset_indPhen_df$phen_temp), scale(temp_pred), method = "kendall")
    Random_cor_salpredict_salphen[i] <- cor(scale(subset_indPhen_df$phen_sal), scale(sal_pred), method = "kendall")
  }
  
  pdf(paste0(path,seed,"_pdf_6zRDA_predict.pdf"), width=10, height=5)
  par(mar=c(1,1))
  plot(nloci, Random_cor_temppredict_tempphen, type="l", ylim=c(-0.1,1), ylab="Cor(RDA prediction, true phenotype)", bty="l" , 
       xlab="Number random loci in RDA", col="darkred", lwd=2, xlim=c(0, 20000), main = paste0(plotmain), cex.main=0.5, las=2)
  points(nloci, Random_cor_salpredict_salphen, type="l", col="cornflowerblue", lwd=2, lty=2)
  legend(10000, 0.1, c("Temp", "Env2"), col=c("darkred", "cornflowerblue"), lty = c(1,2), bty="n")
  
  dev.off()
  
  out_RDA_pred <- data.frame(nloci, Random_cor_temppredict_tempphen, Random_cor_salpredict_salphen)
  
  cor_RDA20000_RDloadings_tempPhen <- Random_cor_temppredict_tempphen[nloci==nmax]
  cor_RDA20000_RDloadings_salPhen <- Random_cor_salpredict_salphen[nloci==nmax]
  
  write.table(out_RDA_pred,paste0(path,seed,"_RDA_predictions.txt"))
  
### AF as a function of environment ####
  
  #str(af_pop)
  #subset_temp_opt
  #subset_sal_opt
  
  sal_levels <- sort(as.numeric(unique(subset_indPhen_df$sal_opt)))
  #sal_levels
  temp_levels <- sort(as.numeric(unique(subset_indPhen_df$temp_opt)))
  #temp_levels
  
  af_sal <- matrix(NA, nrow=length(sal_levels), ncol=nrow(G_full_subset))
  af_temp <- matrix(NA, nrow=length(temp_levels), ncol=nrow(G_full_subset))

  
  # WARNING: this code only works with equal sample sizes for each population
  for (i in 1:ncol(af_pop)){
    af_temp[,i] <- tapply(af_pop[,i], subset_temp_opt, mean)
    af_sal[,i] <- tapply(af_pop[,i], subset_sal_opt, mean)
  }
  colnames(af_temp) <- muts_full$mutID
  v <- tapply(af_pop[,1], subset_temp_opt, mean)
  rownames(af_temp) <- as.numeric(dimnames(v)[[1]])
  
  colnames(af_sal) <- muts_full$mutID
  v <- tapply(af_pop[,1], subset_sal_opt, mean)
  rownames(af_sal) <- as.numeric(dimnames(v)[[1]])

  muts_full$af_cor_temp_pooled <- as.numeric(cor(af_temp, temp_levels))
  muts_full$af_cor_sal_pooled <- as.numeric(cor(af_sal, sal_levels))

  pdf(paste0(path,seed,"_pdf_7_afcors.pdf"), width=5, height=5)
  
  ggplot() + 
    geom_point(data=muts_full[muts_full$causal_temp=="neutral-linked",], aes(x=af_cor_temp, y = af_cor_temp_pooled), color=adjustcolor("goldenrod", 0.1)) +
    geom_point(data=muts_full[muts_full$causal_temp=="neutral",], aes(x=af_cor_temp, y = af_cor_temp_pooled), color=adjustcolor("grey", 0.1)) +
    geom_point(data = muts_full[muts_full$causal_temp=="causal",], 
               aes(x=af_cor_temp, y = af_cor_temp_pooled,  size=Va_temp_prop, color=Va_temp_prop), alpha=0.8, shape=shape_causal) +
    scale_colour_viridis(option=cs, begin = begin_cs, end=end_cs, limits=c(0,1)) +  theme_classic() +  xlab("Cor(af, temp)") + 
    ylab("Cor(af, temp) pooled by temp.") + ggtitle(paste0(plotmain, " mutations: temp")) + ylim(-1,1) + xlim(-1, 1) + geom_abline(data=NULL, intercept=0, slope=1) + labs(color="temp VA prop.", size="temp VA prop.")
    
  ggplot() + 
    geom_point(data=muts_full[muts_full$causal_sal=="neutral-linked",], aes(x=af_cor_sal, y = af_cor_sal_pooled), color=adjustcolor("goldenrod", 0.1)) +
    geom_point(data=muts_full[muts_full$causal_sal=="neutral",], aes(x=af_cor_sal, y = af_cor_sal_pooled), color=adjustcolor("grey", 0.1)) +
    geom_point(data = muts_full[muts_full$causal_sal=="causal",], 
               aes(x=af_cor_sal, y = af_cor_sal_pooled,  size=Va_sal_prop, color=Va_sal_prop), alpha=0.8, shape=shape_causal) +
    scale_colour_viridis(option=cs, begin = begin_cs, end=end_cs, limits=c(0,1)) +  theme_classic() +  xlab("Cor(af, Env2)") + 
    ylab("Cor(af, Env2) pooled by Env2.") + ggtitle(paste0(plotmain, " mutations: Env2")) + ylim(-1,1) + xlim(-1, 1) + geom_abline(data=NULL, intercept=0, slope=1) + labs(color="Env2 VA prop.", size="Env2 VA prop.")

  dev.off()  
  
# Plot phenotype vs. AF clines ####
  
  pdf(paste0(path,seed,"_pdf_8phen-env_af-env.pdf"), width=9, height=9)
  
  par(mfrow=c(2,2), mar=c(4,4,1,1), oma=c(0,0,2,1))

  plot(jitter(subset_indPhen_df$temp_opt), subset_indPhen_df$phen_temp, 
       xlab="Deme Temperature", ylab="Ind. optimum temp.",pch=18,  
       las=1, bty="n", col=adjustcolor("purple",0.1), xlim=c(-1,1), ylim=c(-1.2,1.2))
  abline(lm(subset_indPhen_df$phen_temp~subset_indPhen_df$temp_opt),col="purple", lwd=2)  
  abline(0,1, lty=2)
  
  plot(jitter(subset_indPhen_df$sal_opt), subset_indPhen_df$phen_sal, 
       xlab="Deme Env2", ylab="Ind. optimum Env2.",pch=18,  
       las=1, bty="n", col=adjustcolor("blue",0.1), xlim=c(-1,1), ylim=c(-1.2,1.2))
  abline(lm(subset_indPhen_df$phen_sal~subset_indPhen_df$sal_opt),col="blue", lwd=2)  
  abline(0,1, lty=2)

  mtext(paste0(plotmain), side=3, outer=TRUE, cex=0.5, line=0)

  col <- (mako(10))
  
  #Temp plot
  plot(pop_df$opt1, rep(0, length(pop_df$opt1)), ylim=c(0,1), xlab="Deme Temperature", ylab="Allele frequency", col=rgb(0,0,0,0), bty="n", xlim=c(-1, 1), las=1)
  
  #str(af_temp)
  for (i in which(muts_full$causal_temp=="causal")){
    cor_mut <- cor(temp_levels, af_temp[,i])
    if(abs(cor_mut)<0.15){
      lwd=3
      alpha1 = 0.5
    }else{
      lwd=1
      alpha1 = 0.5
    }
    lines(temp_levels, af_temp[,i], col=adjustcolor(col[ceiling(abs(cor_mut*10))], alpha1), lwd=lwd)
  }  

  plot(pop_df$opt0, rep(0, length(pop_df$opt0)), ylim=c(0,1), xlab="Deme Env2", ylab="Allele frequency", col=rgb(0,0,0,0), bty="n", xlim=c(-1, 1), las=1)  
  
  for (i in which(muts_full$causal_sal=="causal")){
    cor_mut <- cor(sal_levels, af_sal[,i])
    if(abs(cor_mut)<0.15){
      lwd=3
      alpha1 = 0.5
    }else{
      lwd=1
      alpha1 = 0.5
    }
    lines(sal_levels, af_sal[,i], col=adjustcolor(col[ceiling(abs(cor_mut*10))], alpha1), lwd=lwd)
  }  

  par(mfrow=c(1,1))
  plot(pop_df$opt1, rep(0, length(pop_df$opt1)), ylim=c(0,1), xlab="Temperature", ylab="Allele frequency", col=rgb(0,0,0,0), bty="n", xlim=c(-1, 1), main=plotmain)
  
  legend(0,1,seq(0.05,0.95, length.out=10), fill=col, cex=0.7, bty="n", title="Cor (p, env)")  
  legend(-1,0.5, c("linear model", "1:1 line") , lty=c(1,2))
  legend(0.5, 0.5, c("abs(cor) < 0.15", "abs(cor) > 0.15"), lwd=c(3, 0.5))

  dev.off()
  
  
### genotype heatmaps
  
pdf(paste0(path,seed,"_pdf_heatmaps.pdf"), width=5, height=5)

# Temps
order_temp <- order(subset_indPhen_df$temp_opt)
order_sal <- order(subset_indPhen_df$sal_opt)

par(mfrow=c(1,1), mar=c(4,0,6,4))

heatmap(t(G_full_subset[muts_full$causal_temp=="causal",order_temp]), Rowv = NA,  main=paste0(plotmain," Causal alleles Temp"),cexCol = 0.3,  useRaster=FALSE, 
        scale="none",
        #Colv = NA,
        labRow = round(subset_indPhen_df$temp_opt[order_temp],1), xlab="Mutation ID", ylab="South <---------> North")


if(sum(muts_full$causal_sal=="causal")>0){
  heatmap(t(G_full_subset[muts_full$causal_sal=="causal",order_sal]), Rowv = NA,  main=paste0(plotmain,"Causal alleles Env2"),cexCol = 0.3,  useRaster=FALSE,
          scale="none", xlab="Mutation ID", ylab="East <---------> West",
          #Colv = NA,
          labRow = round(subset_indPhen_df$temp_opt[order_temp],1))
}

heatmap(t(G_full_subset[sample(which(!muts_full$causal), 100), order_temp]), Rowv = NA,  main=paste0(plotmain,"Non-causal alleles"),cexCol = 0.3,  useRaster=FALSE,
        scale="none",
        #Colv = NA,
        labRow = round(subset_indPhen_df$temp_opt[order_temp],1), xlab="Mutation ID", ylab="South <---------> North")  


### Visualize genotypes at low salinities ###
subset_indPhen_df2 <- subset_indPhen_df[subset_indPhen_df$sal_opt==-1,]
order_temp2 <- order(subset_indPhen_df2$temp_opt)

G_full_subset2 <- G_full_subset[,subset_indPhen_df$sal_opt==-1]

heatmap(t(G_full_subset2[muts_full$causal,order_temp2]), Rowv = NA,  main=paste0(plotmain," Low Env2. Demes; All causal alleles"),cexCol = 0.3,  useRaster=FALSE, 
        scale="none",
        #Colv = NA,
        labRow = round(subset_indPhen_df2$temp_opt[order_temp2],1), xlab="Mutation ID", ylab="South <---------> North")

#heatmap(t(G_full_subset2[sample(which(!muts_full$causal), 1000),order_temp2]), Rowv = NA,  main=paste0(plotmain," Low Sal. Demes; Non-causal alleles"),cexCol = 0.3,  useRaster=FALSE, 
#        scale="none",
#        #Colv = NA,
#        labRow = round(subset_indPhen_df2$temp_opt[order_temp2],1), xlab="Mutation ID", ylab="South <---------> North")

subset_indPhen_df3 <- subset_indPhen_df[subset_indPhen_df$sal_opt==1,]
order_temp3 <- order(subset_indPhen_df3$temp_opt)

G_full_subset3 <- G_full_subset[,subset_indPhen_df$sal_opt==1]

heatmap(t(G_full_subset3[muts_full$causal,order_temp3]), Rowv = NA,  main="High Env2. Demes; All causal alleles",cexCol = 0.3,  useRaster=FALSE,
        scale="none",
        #Colv = NA,
        labRow = round(subset_indPhen_df3$temp_opt[order_temp3],1), xlab="Mutation ID", ylab="South <---------> North")

# heatmap(t(G_full_subset3[sample(which(!muts_full$causal), 100),order_temp3]), Rowv = NA,  main="High Sal. Demes; Non-causal alleles",cexCol = 0.3,  useRaster=FALSE,
#        scale="none",
#        #Colv = NA,
#        labRow = round(subset_indPhen_df2$temp_opt[order_temp3],1), xlab="Mutation ID", ylab="South <---------> North")

dev.off()

rm(subset_indPhen_df2)
rm(subset_indPhen_df3)
rm(G_full_subset2)
rm(G_full_subset3)


### Write outputs ####

write.table(muts_full, paste0(path,seed,"_muts_full.txt"))
write.table(subset_indPhen_df, paste0(path,seed,"_ind_subset.txt"))


# cor_VA_temp
# cor_VA_sal
# cor_TPR_temp
# cor_FDR_temp
# cor_TPR_sal
# cor_FDR_sal

out_full <- data.frame(seed=as.character(seed), n_samp=n, K,
                       all_corr_phen_temp,
                       subsamp_corr_phen_temp,
                       all_corr_phen_sal,
                       subsamp_corr_phen_sal,
                       num_causal,
                       num_non_causal = sum(!muts_full$causal),
                       num_neut = sum(muts_full$causal_temp=="neutral"), #loci in 2nd half of genome
                       num_causal_temp = sum(muts_full$causal_temp=="causal"),
                       num_causal_sal = sum(muts_full$causal_sal=="causal"),
#                       prop_causal_sig_temp_corr,

                       LEA3.2_lfmm2_Va_temp_prop, 
                       LEA3.2_lfmm2_Va_sal_prop, 
                       LEA3.2_lfmm2_TPR_temp,
                       LEA3.2_lfmm2_TPR_sal,

                      LEA3.2_lfmm2_FDR_allSNPs_temp, 
                      LEA3.2_lfmm2_FDR_allSNPs_sal, 
                      LEA3.2_lfmm2_FDR_neutSNPs_temp,
                      LEA3.2_lfmm2_FDR_neutSNPs_sal,
                      
                      LEA3.2_lfmm2_AUCPR_temp_allSNPs,
                      LEA3.2_lfmm2_AUCPR_temp_neutSNPs,
                      LEA3.2_lfmm2_AUCPR_sal_allSNPs,
                      LEA3.2_lfmm2_AUCPR_sal_neutSNPs,

                       RDA_Va_temp_prop,
                       RDA_Va_sal_prop,
                       RDA_TPR, 
                        RDA_FDR_allSNPs,
                        RDA_FDR_neutSNPs,
                        RDA_AUCPR_allSNPs,
                        RDA_AUCPR_neutSNPs,
                       cor_RDA20000_RDloadings_tempPhen, 
                       cor_RDA20000_RDloadings_salPhen,
 
                      cor_VA_temp_prop,
                      cor_VA_sal_prop,
                      cor_TPR_temp,
                      cor_TPR_sal,
                      
                      cor_FDR_allSNPs_temp,
                      cor_FDR_neutSNPs_temp,
                      cor_FDR_allSNPs_sal,
                      cor_FDR_neutSNPs_sal,
                      
                      cor_AUCPR_temp_allSNPs,
                      cor_AUCPR_temp_neutSNPs,
                      cor_AUCPR_sal_allSNPs,
                      cor_AUCPR_sal_neutSNPs
)


  