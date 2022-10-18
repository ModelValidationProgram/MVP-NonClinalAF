# 
# 
# 
# packages_needed <- c("vcfR", "distances","ggplot2",  "fields", "stringr", "vegan", "robust", "mvtnorm", "viridis", "gridExtra", "devtools")
# 
#  for (i in 1:length(packages_needed)){
#    if(!(packages_needed[i] %in% installed.packages())){install.packages(packages_needed[i])}
#  }
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
                      "PRROC", "qvalue", "OutFLANK", "LEA", "ggExtra")
for (i in 1:length(libraries_needed)){
  library( libraries_needed[i], character.only = TRUE, lib.loc = "/home/lotterhos/R/x86_64-pc-linux-gnu-library/4.0") # for OOD
#  library(libraries_needed[i], character.only = TRUE, lib.loc = "/home/lotterhos/miniconda3/envs/MVP_env_R4.0.3/lib/R/library") # for bash script
}

setwd("/work/lotterhos/MVP-NonClinalAF")

### load params ####

args = commandArgs(trailingOnly=TRUE)

#seed = 1231214 # MAIN GRAPHS
#seed = 1231094
#seed = 1231098
#seed = 1231144
#seed=1232947
#seed = 1231102
#path = "sim_output_20220428/"
#runID=20220428
seed = args[1]
path = args[2]
runID = args[3]

### load data ####

info <- read.table(paste0(path,seed,"_info.txt"), header=TRUE, 
                   colClasses = c("character", 
                                  rep("numeric",15),
                                  "character",
                                  rep("numeric",6)))
info

allsims <- load(paste0("src/0b-final_params-",runID,".RData"))
allsims<- final
thissim <- allsims[grep(seed, allsims$seed),]
(plotmain <- paste(thissim$level, seed, sep="\n"))

vcf_full <- read.vcfR(paste0(path,seed,"_plusneut_MAF01.recode2.vcf.gz"))

vcf_muts <- read.vcfR(paste0(path,seed,"_VCF_causal.vcf.gz"))

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
  
  muts_df$va_temp_full <- muts_df$p*(1-muts_df$p)*muts_df$mutTempEffect^2
  muts_df$va_sal_full <- muts_df$p*(1-muts_df$p)*muts_df$mutSalEffect^2
  #plot(round(vatemp/sum(vatemp),2), abs(muts_df$cor_temp), xlim=c(0,1), ylim=c(0,1))
  #if(sum(vasal)>0){
  #  plot(round(vasal/sum(vasal),2), abs(muts_df$cor_sal), xlim=c(0,1), ylim=c(0,1))
  #}
  va_temp_total <- sum(muts_df$va_temp_full)
  va_sal_total <- sum(muts_df$va_sal_full)
  
  muts_df$va_temp_full_prop <- muts_df$va_temp_full/va_temp_total
  muts_df$va_sal_full_prop <- muts_df$va_sal_full/va_sal_total

  
  


#### plot subpops and migration ####
  #allsim
  mig <- read.table(paste0(path,seed,"_popInfo_m.txt"), header=TRUE)
  
  
  ggtheme <- theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border=element_blank(), 
                                axis.line = element_line(colour="grey20"), axis.title = element_text(colour="grey20"), axis.text = (element_text(colour="grey20")), 
                                legend.title = element_text(colour="grey20"), legend.text = element_text(colour="grey20"))

  
  mig_thick <- rep(NA, length(mig$m))
  mig_thick[mig$m<0.03] <- 0.1
  mig_thick[mig$m>=0.03 & mig$m<0.07] <- 1.5
  mig_thick[mig$m>=0.07 & mig$m<0.15] <- 3
  mig_thick[mig$m>=0.15 & mig$m<0.3] <- 5
  mig_thick[mig$m>=0.3] <- 7
  
  
  par(mfrow=c(1,1), mar=c(4,4,3,1))
  pdf(paste0(path,seed,"_pdf_1pop_mig.pdf"), width=8, height=8)  
  plot(pop_df$x, pop_df$y, col = rgb(0,0,0,0), bty="l", xlab="x", ylab="y", main=plotmain, cex.main=0.5)
  for (i in 1:nrow(mig)){
    start_x <- pop_df$x[which(pop_df$subpopID==mig$from[i])]
    start_y <- pop_df$y[which(pop_df$subpopID==mig$from[i])]
    end_x <- pop_df$x[which(pop_df$subpopID==mig$to[i])]
    end_y <- pop_df$y[which(pop_df$subpopID==mig$to[i])]
    adj = 0.01  
    if (end_x < start_x){end_x <- end_x +0.3; start_x <- start_x - 0.5}
    if (end_x > start_x){end_x <- end_x -0.3; start_x <- start_x + 0.5}
    if (end_y < start_y){end_y <- end_y +0.3; start_y <- start_y - 0.5}
    if (end_y > start_y){end_y <- end_y -0.3; start_y <- start_y + 0.5}
    arrows(start_x,start_y,end_x, end_y, col="cornflowerblue", lwd=mig_thick[i], length=0.05)
  }
  text(pop_df$x, pop_df$y, pop_df$N, cex=1)
  dev.off()
  
  pdf(paste0(path,seed,"_pdf_1pop.pdf"), width=5, height=5)
  ggplot(pop_df) + ggtheme + geom_point(aes(x=x, y=y,size=opt0), color="grey20") + geom_point(aes(x=x, y=y, color=opt1), size=2.5) + scale_colour_gradient2(high=rgb(1,0.4,0.2), low="cornflowerblue", mid=rgb(0.8,0.8,0.7), name="Temp") + geom_text(aes(x=x, y=y+0.3,label=subpopID)) + labs(size="Env2") + ggtitle(plotmain) + theme(plot.title = element_text(size = 7))
  ggplot(pop_df) + ggtheme + geom_point(aes(x=x, y=y,size=opt0), color="grey20") + geom_point(aes(x=x, y=y, color=opt1), size=2.5) + scale_colour_gradient2(high=rgb(1,0.4,0.2), low="cornflowerblue", mid=rgb(0.8,0.8,0.7), name="Temp") +  labs(size="Env2") + ggtitle(plotmain) + theme(plot.title = element_text(size = 7))
  
  #pdf(paste0(seed,"_pdf_1apop.pdf"), width=5, height=4)
  #f <- ggplot(pop_df) + geom_point(aes(x=x, y=y,size=opt0)) + geom_point(aes(x=x, y=y, color=opt1)) + theme(plot.background = element_rect(fill = "transparent",colour = NA)) +  scale_colour_gradient2(high=rgb(1,0.4,0.2), low="cornflowerblue", mid=rgb(0.8,0.8,0.7), name="Temp") +  labs(size="Env2") + ggtitle(plotmain) + theme(plot.title = element_text(size = 7))
  #ggsave(paste0(seed,"_pdf_1apop.png"), f, bg="transparent")
  
  plot(LA_df$gen, LA_df$cor_temp_ind, type="l", col="orange", lwd=3, main=plotmain, 
       ylab="cor (Env, phenotype)", xlab="Generation", ylim=c(0,1), bty="n")
  
  if (sum(LA_df$cor_sal_ind, na.rm=TRUE)>0){
    points(LA_df$gen, LA_df$cor_sal_ind, col="cornflowerblue", type="l", lwd=3, lty=2)
  }
  legend(0,1, c("Temp", "Env2"), lwd=3, col=c("orange", "cornflowerblue"), lty=c(1,2), bty="n")
  
  plot(LA_df$gen, LA_df$local_adapt, type="l", col="blue", main=plotmain, 
       ylab="Amount of local adaptation", xlab="Generation", ylim=c(0,1.2), bty="n")
  
  plot(LA_df$gen, LA_df$mean_phen1, type="l", col="blue", main=plotmain, 
       ylab="Mean metapopulation phenotype", xlab="Generation", ylim=c(-1,1))
    if (sum(LA_df$cor_sal_ind, na.rm=TRUE)>0){
      points(LA_df$gen, LA_df$mean_phen0, col="orange", type="l", lwd=3, lty=2)
    }
  
  # Find high fitness individuals ####
  npops <- length(levels(factor(indPhen_df$subpop)))
  n = 10 # number of individual per pop
  n_per_pop = n # for output
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
  
  boxplot(indPhen_df$fitness[indPhen_df$subset]~indPhen_df$N[indPhen_df$subset], xlab="N", ylab="fitness", ylim=c(0,1))
  plot(indPhen_df$fitness[indPhen_df$subset]~indPhen_df$N[indPhen_df$subset], xlab="N", ylab="fitness", ylim=c(0,1))
  # looks good to me
  
  hist_fitness <- tapply(indPhen_df$fitness[indPhen_df$subset],indPhen_df$subpopID[indPhen_df$subset], mean)
  sd_fitness <- tapply(indPhen_df$fitness[indPhen_df$subset],indPhen_df$subpopID[indPhen_df$subset], sd)
  
  sd_fitness_among_inds <- sd(indPhen_df$fitness[indPhen_df$subset])
  sd_fitness_among_pops <- sd(hist_fitness) 
  
  hist(as.numeric(hist_fitness), breaks=seq(0,1, 0.01), main=plotmain, cex.main=0.5)
  hist(as.numeric(sd_fitness), breaks=seq(0,1, 0.01), main=plotmain, cex.main=0.5)
  
  dev.off()
  ## end plot high fitness individuals ####
  
  if (var(indPhen_df$phen_sal)>0){
    all_corr_phen_sal <- cor(indPhen_df$phen_sal,indPhen_df$sal_opt)
    subsamp_corr_phen_sal <- cor(indPhen_df$phen_sal[which(indPhen_df$subset)], indPhen_df$sal_opt[which(indPhen_df$subset)], method="kendall")
  }else{
    all_corr_phen_sal <- NA
    subsamp_corr_phen_sal <- NA
  }
  
  all_corr_phen_temp <- cor(indPhen_df$phen_temp,indPhen_df$temp_opt)
  
  subsamp_corr_phen_temp <- cor(indPhen_df$phen_temp[which(indPhen_df$subset)],indPhen_df$temp_opt[which(indPhen_df$subset)], method="kendall")

### show cool patterns of mutation ####
  
  pdf(paste0(path,seed,"_pdf_2muts.pdf"), width=5, height=5)
  
  plot(muts_df$mutSalEffect, muts_df$mutTempEffect, ylim=c(-0.7, 0.7), xlim=c(-0.7,0.7), xlab="mutation effect on Env2 phenotype", 
       ylab="mutation effect on temperature phenotype", bty="l", main=plotmain, cex.main=0.5)
  abline(h=0, col="grey")
  abline(v=0, col="grey")
  
  geno_full <- vcf_full@gt[,-1] 
  dim(geno_full)
  position_full <- getPOS(vcf_full)
  rown <- vcf_full@fix[,"ALT"] # mutation ID
  # rown = 1 is a neutral mutation
  causal_mut_locs <- which(rown %in% as.character(muts_df$mutID)) # causal mutations
  #position_full[causal_mut_locs] # positions of causal mutations
  
  
  head(geno_full[,1:5])
  
  
  G_full <- matrix(NA, nrow = nrow(geno_full), ncol = ncol(geno_full))
  G_full[geno_full %in% c("0/0", "0|0")] <- 0
  G_full[geno_full  %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
  G_full[geno_full %in% c("1/1", "1|1")] <- 2

  
  
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
  
  # Assign loci to LG groups
  numLGs <- ceiling(max(as.numeric(vcf_full@fix[,"POS"]))/50000)
  for (LG in 1:numLGs){
    sites_low <- (LG-1)*50000+1
    sites_high <- (LG)*50000
    pos <- as.numeric(vcf_full@fix[,"POS"])
    whichlocs <- pos >= sites_low & pos <= sites_high
    vcf_full@fix[whichlocs,"CHROM"] <- as.numeric(LG)
  }
  table(vcf_full@fix[,"CHROM"])
  
  # Assign unique names to each locus
  head(vcf_full@fix)
  vcf_full@fix[,"INFO"] <- paste0(vcf_full@fix[,"CHROM"],"-",vcf_full@fix[,"POS"])
  if (sum(duplicated(vcf_full@fix[,"POS"]))>0){
    dups <- which(duplicated(vcf_full@fix[,"POS"]))
    vcf_full@fix[dups,"INFO"] <- paste0(vcf_full@fix[dups,"INFO"],"_2")
  }
  b <- 3
  while (sum(duplicated(vcf_full@fix[,"INFO"]))>0){
    dups <- which(duplicated(vcf_full@fix[,"INFO"]))
    vcf_full@fix[dups,"INFO"] <- gsub(paste0("_",b-1),
                                      paste0("_",b),
                                      x = vcf_full@fix[dups,"INFO"])
    b <- b+1
    if (b==10){print("problem in naming loci"); break}
  }
  
  # Check duplicated locus names
  if(sum(duplicated(vcf_full@fix[dups,"INFO"]))>0){print("Error:duplicate locus names");break}
  
  head(vcf_full@fix)
  
  rownames(G_full_subset)=vcf_full@fix[,"INFO"]
  
  colnames(G_full_subset)=indPhen_df$indID[which(indPhen_df$subset)]
  head(G_full_subset[,1:5])
  plot(colnames(G_full_subset), indPhen_df$indID[which(indPhen_df$subset)])
  
  
  a_freq_subset <- as.numeric(rowSums(G_full_subset)/(2*ncol(G_full_subset)))
  # should have 1000 individuals
  #hist(a_freq_subset)
  
  muts_full0 <- data.frame(seed = as.character(seed),
                           VCFrow = 1:nrow(vcf_full@gt),
                           mutID = vcf_full@fix[,"ALT"],
                           LG = as.integer(vcf_full@fix[,"CHROM"]),
                           pos_pyslim = as.integer(vcf_full@fix[,"POS"]),
                           mutname = vcf_full@fix[,"INFO"],
                           a_freq_full = a_freq_full,
                           a_freq_subset = a_freq_subset)
  
  if(!identical(rownames(G_full_subset), muts_full0$mutname)){Print("Error mut names out of order");break}
  
  head(muts_full0)
  dim(muts_full0)
  muts_full0$mutID <- as.character(muts_full0$mutID) # all loci, but only causal loci have unique ID
  muts_df$mutID <- as.character(muts_df$mutID) #only causal loci, unique ID from SliM
  
  muts_full <- merge(muts_full0, muts_df, by=c("mutID", "seed"), all.x=TRUE)
  muts_full <- muts_full[order(muts_full$VCFrow),]
  head(muts_full)
  
  if(!identical(rownames(G_full_subset), muts_full$mutname)){Print("Error mut names out of order");break}
  
  
  # Are there any missing values for LG?
  if(sum(is.na(muts_full$LG))>0){Print("Error: Missing LGs"); break}
  
  # make sure mutations in correct order
  #plot(muts_full$pos_pyslim)
  rm(muts_full0)
  rm(muts_df)
  dim(muts_full)
  
  muts_full$mutTempEffect <- as.numeric(muts_full$mutTempEffect)
  muts_full$mutSalEffect <- as.numeric(muts_full$mutSalEffect)
  muts_full$causal <- FALSE
  muts_full$causal[muts_full$mutID!=1] <- TRUE
  
  
  # FILTER OUT ALLELS WITH MAF < 0.01 AFTER SAMPLING ####
   if(!identical(as.character(muts_full$mutname), as.character(rownames(G_full_subset)))){print("Error 1a: mutations not lined up");break()}
    
    #muts_full$isMAF01 <- FALSE
    rm_loci <- which(muts_full$a_freq_subset < 0.01 | muts_full$a_freq_subset > 0.99)
    #muts_full$isMAF01[rm_loci] <- rep(TRUE, length=length(rm_loci))
  
    num_multiallelic <- sum((!complete.cases(G_full_subset)))
    rm_loci <- sort(c(rm_loci, which(!complete.cases(G_full_subset))))
        
    numCausalLowMAFsample <- sum(muts_full$a_freq_subset < 0.01 & muts_full$causal)
    
    # Filter out MAF loci from G_full_subset
    
    if(length(rm_loci)>0){
      G_full_subset <- G_full_subset[-rm_loci,]
      muts_full <- muts_full[-rm_loci,]
    }
    if(!identical(as.character(muts_full$mutname), as.character(rownames(G_full_subset)))){print("Error 1a2: mutations not lined up");break()}
    
    dim(G_full_subset)
    
    head(G_full_subset[,1:10])
    
    hist(muts_full$a_freq_subset, breaks=seq(0,1,0.01))
  plot(muts_full$a_freq_full, muts_full$a_freq_subset)
  hist(muts_full$mutTempEffect, xlim=c(-2,2), xlab="Mutation temp effect", breaks=seq(-5,5,0.01), main=plotmain, cex.main=0.5)
  hist(muts_full$mutSalEffect, xlim=c(-2,2), xlab="Mutation Env2 effect", breaks=seq(-5,5,0.01), main=plotmain, cex.main=0.5)
  
  # Sanity check that muts_full mutations line up with G_full_subset
  if(nrow(muts_full)!=nrow(G_full_subset)){print("Error 1: mutation data frames have different number of rows");break()}

  num_causal_postfilter <- sum(muts_full$causal)
  num_causal_prefilter <- nrow(vcf_muts@fix)
  num_neut_postfilter <- sum(!muts_full$causal)
  num_neut_prefilter <- sum((vcf_full@fix[,"ALT"]=="1"))
  
  dim(G_full_subset)
  
  subset_indPhen_df <- indPhen_df[indPhen_df$subset,]
  if(!identical(subset_indPhen_df$indID, sort(subset_indPhen_df$indID))){print("Error individuals in wrong order"); break()}

  str(subset_indPhen_df)
  
  head(pop_df)
  n_ind <- table(subset_indPhen_df$subpop)
  

  subset_temp_opt <- tapply(subset_indPhen_df$temp_opt,
                            subset_indPhen_df$subpop, mean)
  
  subset_sal_opt <- tapply(subset_indPhen_df$sal_opt,
                           subset_indPhen_df$subpop, mean)
  
  subpop_subset <- as.numeric(rownames(subset_sal_opt))
  
  if(!identical(as.character(muts_full$mutname), as.character(rownames(G_full_subset)))){print("Error 1b: mutations not lined up");break()}
  
  
  af_pop <- matrix(NA, nrow=length(subpop_subset), ncol=nrow(G_full_subset))
  
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
  #plotsamp <- sample(1:nrow(muts_full), 2000, replace=FALSE)
  plot(muts_full$cor_temp, muts_full$af_cor_temp, ylim=c(-1,1), xlim=c(-1,1), xlab="SliM cor(p, Temp)", ylab="Subsample cor(p, Temp)", main=plotmain)
  abline(0,1)
  
  plot(muts_full$cor_sal, muts_full$af_cor_sal, ylim=c(-1,1), xlim=c(-1,1), xlab="SliM cor(p, Env2)", ylab="Subsample cor(p, Env2)", main=plotmain)
  abline(0,1)
  
  plot(muts_full$cor_temp, muts_full$cor_sal, main="SliM Correlation (AF, ENV) all samples", xlim=c(-1,1), ylim=c(-1,1))
  plot(muts_full$af_cor_temp, muts_full$af_cor_sal, main="Subsample Correlation (AF, ENV) all samples", xlim=c(-1,1), ylim=c(-1,1))
  
  if(!identical(as.character(muts_full$mutname), as.character(rownames(G_full_subset)))){print("Error 1c: mutations not lined up");break()}
  
  
  
  muts_full$causal_temp <- "neutral-linked"
  muts_full$causal_temp[muts_full$causal & muts_full$mutTempEffect != 0] <- "causal"
  muts_full$causal_temp[muts_full$pos_pyslim > 500000] <- "neutral"
  

  muts_full$causal_sal <- "neutral-linked"
  muts_full$causal_sal[muts_full$causal & muts_full$mutSalEffect != 0 & info$Ntraits == 2] <- "causal"
  muts_full$causal_sal[muts_full$pos_pyslim > 500000] <- "neutral"
  

  
  muts_full$pos_pyslim <- as.integer(muts_full$pos_pyslim)
  
  muts_full$colors <- NA
  muts_full$colors[muts_full$LG %in% seq(1,9, by=2)] <- adjustcolor("goldenrod1", 0.1)
  muts_full$colors[muts_full$LG %in% seq(2,10, by=2)] <- adjustcolor("goldenrod3", 0.1)
  muts_full$colors[muts_full$LG %in% seq(11,19, by=2)] <- adjustcolor("grey50", 0.1)
  muts_full$colors[muts_full$LG %in% seq(12,20, by=2)] <- adjustcolor("grey70", 0.1)
  muts_full$colors <- as.character(muts_full$colors)
  
  muts_full$Va_temp <- muts_full$a_freq_subset*(1-muts_full$a_freq_subset)*muts_full$mutTempEffect^2
  muts_full$Va_temp[is.na(muts_full$Va_temp)] <- 0
  muts_full$Va_temp_prop <- muts_full$Va_temp/sum(muts_full$Va_temp, na.rm=TRUE)
  
  hist(muts_full$Va_temp_prop[muts_full$causal_temp=="causal"], xlab="VA temp prop", main=plotmain, xlim=c(0,1), breaks=seq(0,1,0.01))
  muts_full$Va_temp_prop[is.na(muts_full$Va_temp_prop)] <- 0
  
  muts_full$Va_sal <- muts_full$a_freq_subset*(1-muts_full$a_freq_subset)*muts_full$mutSalEffect^2
  muts_full$Va_sal[is.na(muts_full$Va_sal)] <- 0
  muts_full$Va_sal_prop <- muts_full$Va_sal/sum(muts_full$Va_sal, na.rm=TRUE)
  muts_full$Va_sal_prop[is.na(muts_full$Va_sal_prop)] <- 0
  hist(muts_full$Va_sal_prop[muts_full$causal_sal=="causal"], xlab= "VA Env2 prop.", main= plotmain, xlim=c(0,1), breaks=seq(0,1,0.01))
  
  Va_temp_sample <- sum(muts_full$Va_temp)
  Va_sal_sample <- sum(muts_full$Va_sal)
  
  plot(muts_full$va_temp_full_prop, muts_full$Va_temp_prop)
  abline(0,1)
  
  if(sum(va_sal_total>0)){
    plot(muts_full$va_sal_full_prop, muts_full$Va_sal_prop)
    abline(0,1)
  }

  if(!identical(as.character(muts_full$mutname), as.character(rownames(G_full_subset)))){print("Error 1d: mutations not lined up");break()}
  
  plot(muts_full$Va_temp_prop[muts_full$causal==TRUE],
       muts_full$Va_sal_prop[muts_full$causal==TRUE],
       xlim=c(-0.01,1), ylim=c(-0.01,1), bty="l", pch=19, 
       xlab="Subset VA temp prop.",
       ylab="Subset VA Env2 prop.", main=plotmain)
  
  
  muts_full$causal_temp <- factor(muts_full$causal_temp, levels=c("causal", "neutral", "neutral-linked"),
                                  ordered=TRUE)
  muts_full$causal_sal <- factor(muts_full$causal_sal, levels=c("causal", "neutral", "neutral-linked"),
                                 ordered=TRUE)
  

  ggplot(muts_full, aes(af_cor_temp, fill=causal_temp)) + ggtheme +
    scale_fill_manual(values=c("blue", adjustcolor("grey80", alpha.f=0.4), 
                                adjustcolor("goldenrod1", alpha.f = 0.1)),
                       drop=FALSE) + 
    theme(plot.title=element_text(size=8)) +
    geom_density(alpha=0.5) + xlab("Cor(p, temp)") + xlim(-1,1) + ggtitle(plotmain) 
  
  ggplot(muts_full, aes(af_cor_sal, fill=causal_sal)) + ggtheme +
    scale_fill_manual(values=c("blue", adjustcolor("grey80", alpha.f=0.4), 
                                adjustcolor("goldenrod1", alpha.f = 0.1)),
                       drop=FALSE) +
    theme(plot.title=element_text(size=8)) +
    geom_density(alpha=0.5) + xlab("Cor(p, Env2)") + xlim(-1,1) + 
    ggtitle(plotmain)
  
  
  p1<- ggplot(aes(x=Va_temp_prop,y=abs(af_cor_temp),  fill=causal_temp,
                  shape=causal_temp, colour=causal_temp), data=muts_full) + 
    geom_point() + ggtheme +
    ggtitle(plotmain) + 
    theme(plot.title=element_text(size=8)) +
    scale_x_continuous(limits=c(0,1)) +
    scale_y_continuous(limits=c(0,1)) +
    scale_shape_manual(values=c(3, 0, 1)) +
    scale_color_manual(values=c("blue", adjustcolor("grey60", alpha.f=0.8), 
                                adjustcolor("goldenrod1", alpha.f = 1)),
                       drop=FALSE)
  p1a<-ggExtra::ggMarginal(p1, type="density",groupColour = TRUE, groupFill = TRUE, bw=0.05, size=3)

  

  p<- ggplot(aes(x=Va_sal_prop,y=abs(af_cor_sal),  fill=causal_sal,
                 shape=causal_sal, colour=causal_sal), data=muts_full) + 
    geom_point() + ggtheme +
    ggtitle(plotmain) + 
    theme(plot.title=element_text(size=8)) +
    scale_x_continuous(limits=c(0,1)) +
    scale_y_continuous(limits=c(0,1)) +
    scale_shape_manual(values=c(3, 0, 1), drop=FALSE) +
    scale_color_manual(values=c("blue", adjustcolor("grey60", alpha.f=0.8), 
                                adjustcolor("goldenrod1", alpha.f = 1)),
                       drop=FALSE)
  pqb <- ggExtra::ggMarginal(p, type="density",groupColour = TRUE, groupFill = TRUE, bw=0.05, size=3)
  
  grid.arrange(p1a, pqb, nrow=2)
  
  dev.off()
  
  rm(geno_full)
  rm(indPhen_df)
  rm(G_full)

  ### end show cool patterns of mutations ###

  ## Write G_full_subset to file ####
  if(!identical(as.character(muts_full$mutname), as.character(rownames(G_full_subset)))){print("Error 1e: mutations not lined up");break()}
  
  write.table(G_full_subset, paste0(path,seed,"_Rout_Gmat_sample.txt"), row.names=TRUE, col.names=TRUE)
    print("wrote G matrix to file")
    
    dim(G_full_subset)
    dim(muts_full)
  
  ### Correlation stats ####
  Bonf_alpha <- (0.05/(nrow(muts_full)))
  
  muts_full$cor_temp_sig <- muts_full$af_cor_temp_P < Bonf_alpha
  muts_full$cor_sal_sig <- muts_full$af_cor_sal_P < Bonf_alpha
  
  cor_af_temp_noutliers <- sum(muts_full$cor_temp_sig)
  cor_af_sal_noutliers <- sum(muts_full$cor_sal_sig)
  nSNPs <- nrow(muts_full)
  
  num_causal_sig_temp_corr <- sum(muts_full$af_cor_temp_P[muts_full$causal_temp=="causal"]<Bonf_alpha)# number of causal loci that have significant Kendall's correlations with temperature after Bonferroni correction
  
  num_causal_sig_sal_corr<- sum(muts_full$af_cor_sal_P[muts_full$causal_sal=="causal"]<Bonf_alpha)# number of causal loci that have significant Kendall's correlations with salinity after Bonferroni correction
  
  num_notCausal_sig_temp_corr <- sum(muts_full$af_cor_temp_P[!muts_full$causal_temp=="causal"]<Bonf_alpha)# number of non-causal (neutral and neutral-linked) loci that have significant  Kendall's correlations with temperature after Bonferroni correction
  num_notCausal_sig_sal_corr <- sum(muts_full$af_cor_sal_P[!muts_full$causal_sal=="causal"]<Bonf_alpha)# number of non-causal (neutral and neutral-linked) loci that have significant Kendall's correlations with salinity after Bonferroni correction
  
  num_neut_sig_temp_corr <- sum(muts_full$af_cor_temp_P[muts_full$causal_temp=="neutral"]<Bonf_alpha)# number of neutral (unlinked to causal) loci that have significant  Kendall's correlations with temperature after Bonferroni correction
  num_neut_sig_sal_corr <- sum(muts_full$af_cor_sal_P[muts_full$causal_sal=="neutral"]<Bonf_alpha)# number of neutral (unlinked to causal) loci that have significant Kendall's correlations with salinity after Bonferroni correction
  
  cor_VA_temp_prop <- sum(muts_full$Va_temp_prop[which(muts_full$causal_temp=="causal" & (muts_full$af_cor_temp_P< Bonf_alpha) )]) 
  cor_VA_sal_prop <- sum(muts_full$Va_sal_prop[which(muts_full$causal_sal=="causal" & (muts_full$af_cor_sal_P< Bonf_alpha) )]) 
  
  cor_TPR_temp <- num_causal_sig_temp_corr/sum(muts_full$causal_temp=="causal")
  cor_FPR_temp_neutSNPs <- (num_neut_sig_temp_corr)/sum(muts_full$causal_temp=="neutral")
  cor_FDR_allSNPs_temp <- num_notCausal_sig_temp_corr/sum(muts_full$af_cor_temp_P<Bonf_alpha) # all neutral / all outliers
  cor_FDR_neutSNPs_temp <- num_neut_sig_temp_corr/(num_neut_sig_temp_corr+num_causal_sig_temp_corr) # definitely neutral/ relevant outliers
  
  cor_TPR_sal <- num_causal_sig_sal_corr/sum(muts_full$causal_sal=="causal")
  cor_FPR_sal_neutSNPs <- (num_neut_sig_sal_corr)/sum(muts_full$causal_sal=="neutral")
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
  
  
  if(!identical(as.character(muts_full$mutname), as.character(rownames(G_full_subset)))){print("Error 1f: mutations not lined up");break()}
  
  ### Correlation Manhattan ####
  
  pdf(paste0(path,seed,"_pdf_8manhattan_cor.pdf"), width=15, height=10)
  
  cs <- "mako"
  begin_cs <- 0.9
  end_cs <- 0.2
  shape_causal <- 17
  outlier_color <- adjustcolor("brown1", 0.8)
  ymax=max(c(-log10(muts_full$af_cor_sal_P), -log10(muts_full$af_cor_temp_P), 10))
  
  a <- ggplot() + ggtheme + 
    geom_point(data=muts_full, aes(x=pos_pyslim, y = -log10(af_cor_temp_P)), color=muts_full$colors) +
    geom_point(data = muts_full[muts_full$cor_temp_sig,], aes(x=pos_pyslim, y =-log10(af_cor_temp_P)),pch=23, col=outlier_color, size=3, alpha=0.5) + 
    ylim(0, ymax)  + ylab("-log10(P) Cor(af, temp)") +
    geom_point(data = muts_full[muts_full$causal_temp=="causal",], aes(x=pos_pyslim, y = -log10(af_cor_temp_P), 
                                                                       color=Va_temp_prop, size=Va_temp_prop), shape=shape_causal, alpha=0.7) + 
    scale_colour_viridis(option=cs, begin = begin_cs, end=end_cs, limits=c(0,1)) + 
    ggtitle(paste0(plotmain," temp")) + xlab("position") + labs(color="temp VA prop.", size="temp VA prop.")
  
  b <- ggplot() + ggtheme +
    geom_point(data=muts_full, aes(x=pos_pyslim, y = -log10(af_cor_sal_P)), color=muts_full$colors) +
    geom_point(data = muts_full[muts_full$cor_sal_sig,], aes(x=pos_pyslim, y =-log10(af_cor_sal_P)),pch=23, col=outlier_color, size=3, alpha=0.5) + 
    ylim(0, ymax)  + ylab("-log10(P) Cor(af, Env2)") +
    geom_point(data = muts_full[muts_full$causal_sal=="causal",], aes(x=pos_pyslim, y = -log10(af_cor_sal_P), 
                                                                      color=Va_sal_prop, size=Va_sal_prop), shape=shape_causal, alpha=0.7) + 
    scale_colour_viridis(option=cs, begin = begin_cs, end=end_cs, limits=c(0,1)) + 
    ggtitle(paste0(plotmain," temp")) + xlab("position") + labs(color="Env2 VA prop.", size="Env2 VA prop.")
  
  grid.arrange(a, b, nrow=2)
  
  dev.off()
  
### PCA ####
  lfmmfile <- paste0(path, seed, "_genotypes.lfmm")
  write.lfmm(t(G_full_subset), lfmmfile)
  
  pc = pca(lfmmfile, 30, scale = TRUE)
  saveRDS(pc,paste0(path, seed, "_pca.RDS"))
  print("calculated pca")
  subset_indPhen_df$PC1 <- pc$projections[,1]
  subset_indPhen_df$PC2 <- pc$projections[,2]
  subset_indPhen_df$PC3 <- pc$projections[,3]
  tw = tracy.widom(pc)
  plot(tw$percentage)
  a <- tw$percentage[1:20]
  b <- tw$percentage[2:21]
  K1 <- max(which(a > b*1.5))
  K2 <- max(which(a > b*1.4))
  K3 <- max(which(a > b*1.3)) # a little safety net if the first criteria doesn't work
  if (!is.infinite(K1)){K=K1}
  if (is.infinite(K1)){K=K3}
  if (is.infinite(K)){print("Error K is not definite"); break}
  print(c("K=",K))
  print(c("K1=", K1, "K2=", K2, "K3=", K3))
  
  pdf(paste0(path,seed,"_pdf_3pca.pdf"), width=9, height=8)
  
  propvarpc <- pc$sdev[1:15]/sum(pc$sdev)
  plot(propvarpc, bty="l", ylab="Prop Var of PC axis", main=paste0(plotmain, "; K=", K), cex.main=0.5)
  
  ggplot(subset_indPhen_df) + ggtheme + geom_point(aes(x=PC1, y=PC2,size=sal_opt), color="grey20") + 
    geom_point(aes(x=PC1, y=PC2, color=temp_opt), size=2.5) +  
    scale_colour_gradient2(high=rgb(1,0.4,0.2), low="cornflowerblue", mid=rgb(0.8,0.8,0.7), name="Temp") + 
    labs(size="Env2") + ggtitle(paste0(plotmain, "; K=", K)) + 
    xlab(paste0("PC1 (", round(propvarpc[1]*100,1),"%)")) +
    ylab(paste0("PC2 (", round(propvarpc[2]*100,1),"%)"))
    
  
  ggplot(subset_indPhen_df) + ggtheme + geom_point(aes(x=PC1, y=PC3,size=sal_opt), color="grey20") + geom_point(aes(x=PC1, y=PC3, color=temp_opt), size=2.5) + scale_colour_gradient2(high=rgb(1,0.4,0.2), low="cornflowerblue", mid=rgb(0.8,0.8,0.7), name="Temp") + labs(size="Env2") + ggtitle(paste0(plotmain,  "; K=", K)) + 
    xlab(paste0("PC1 (", round(propvarpc[1]*100,1),"%)")) +
    ylab(paste0("PC3 (", round(propvarpc[3]*100,1),"%)"))
  
  
  dev.off()
  

# FST Outflank ####
  fst <- MakeDiploidFSTMat(t(G_full_subset), locusNames = muts_full$mutname, 
                           popNames = subset_indPhen_df$subpopID)
  
  head(fst)
  
  if(!identical(as.character(fst$locusNames), as.character(rownames(muts_full$mutname)))){print("Error 1h: mutations not lined up");break()}
  
  if(!identical(as.character(muts_full$mutname), as.character(rownames(G_full_subset)))){print("Error 1g: mutations not lined up");break()}
  
  
  muts_full$He_outflank <- fst$He
  muts_full$Fst_outflank <- fst$FST
  meanFst <- mean(fst$T1)/mean(fst$T2)
  
  # FST manhattan ###
  
  pdf(paste0(path,seed,"_pdf_8manhattan_fst.pdf"), width=15, height=5)
  
  cs <- "mako"
  begin_cs <- 0.9
  end_cs <- 0.2
  shape_causal <- 17
  outlier_color <- adjustcolor("brown1", 0.8)
  ymax=max(c(-log10(muts_full$af_cor_sal_P), -log10(muts_full$af_cor_temp_P), 10))
  
  ggplot() + ggtheme + 
    geom_point(data=muts_full, aes(x=pos_pyslim, y = Fst_outflank), color=muts_full$colors) +
    #geom_point(data = muts_full[muts_full$cor_temp_sig,], aes(x=pos_pyslim, y =-log10(af_cor_temp_P)),pch=23, col=outlier_color, size=3, alpha=0.5) + 
    ylim(0, 1)  + ylab("FST") +
    geom_point(data = muts_full[muts_full$causal_temp=="causal",], aes(x=pos_pyslim, y = Fst_outflank, 
                                                                       color=Va_temp_prop, size=Va_temp_prop), shape=shape_causal, alpha=0.7) + 
    scale_colour_viridis(option=cs, begin = begin_cs, end=end_cs, limits=c(0,1)) + 
    ggtitle(paste0(plotmain," temp")) + xlab("position") + labs(color="temp VA prop.", size="temp VA prop.") +
    geom_hline(aes(yintercept=meanFst), color=adjustcolor("darkred",0.5))
  
  dev.off()
  
  print("FST finished")
  
# LFMM ####
  
  mod_temp <- lfmm2(input = t(G_full_subset), env = subset_indPhen_df$temp_opt, K = K)
  
  saveRDS(mod_temp , paste0(path,seed,"_lfmm2_temp.RDS"))
  
  pv_temp <- lfmm2.test(object = mod_temp, input = t(G_full_subset),
                        env = subset_indPhen_df$temp_opt,
                        linear = TRUE)
  
  muts_full$LEA3.2_lfmm2_mlog10P_tempenv <- -log10(as.numeric(pv_temp$pvalues))
  muts_full$LEA3.2_lfmm2_mlog10P_tempenv_sig <- qvalue(as.numeric(pv_temp$pvalues))$qvalue < 0.05
  
  LEA3.2_lfmm2_mlog10P_tempenv_noutliers <- sum(muts_full$LEA3.2_lfmm2_mlog10P_tempenv_sig)
  
  mod_sal <- lfmm2(input = t(G_full_subset), env = subset_indPhen_df$sal_opt, K = K)
  
  saveRDS(mod_sal , paste0(path,seed,"_lfmm2_sal.RDS"))
  
  pv_sal <- lfmm2.test(object = mod_sal, input = t(G_full_subset),
                       env = subset_indPhen_df$sal_opt,
                       linear = TRUE)
  
  
  muts_full$LEA3.2_lfmm2_mlog10P_salenv <- -log10(as.numeric(pv_sal$pvalues))
  muts_full$LEA3.2_lfmm2_mlog10P_salenv_sig <- qvalue(as.numeric(pv_sal$pvalues))$qvalue < 0.05
  
  LEA3.2_lfmm2_mlog10P_salenv_noutliers <- sum(muts_full$LEA3.2_lfmm2_mlog10P_salenv_sig)
  
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
  
  LEA3.2_lfmm2_FPR_neutSNPs_temp <-LEA3.2_lfmm2_num_neut_sig_temp/(sum(muts_full$causal_temp=="neutral"))
  LEA3.2_lfmm2_FPR_neutSNPs_sal <- LEA3.2_lfmm2_num_neut_sig_sal/(sum(muts_full$causal_sal=="neutral"))
  # neutral / relevant outliers
  
  LEA3.2_lfmm2_AUCPR_temp_allSNPs<- (pr.curve(scores.class0 =muts_full$LEA3.2_lfmm2_mlog10P_salenv[muts_full$causal_temp=="causal"], scores.class1 = muts_full$LEA3.2_lfmm2_mlog10P_salenv[!muts_full$causal_temp=="causal"]))$auc.integral
  LEA3.2_lfmm2_AUCPR_temp_neutSNPs<- (pr.curve(scores.class0 =muts_full$LEA3.2_lfmm2_mlog10P_salenv[muts_full$causal_temp=="causal"], scores.class1 = muts_full$LEA3.2_lfmm2_mlog10P_salenv[muts_full$causal_temp=="neutral"]))$auc.integral
  LEA3.2_lfmm2_AUCPR_sal_allSNPs<- (pr.curve(scores.class0 =muts_full$LEA3.2_lfmm2_mlog10P_salenv[muts_full$causal_sal=="causal"], scores.class1 = muts_full$LEA3.2_lfmm2_mlog10P_salenv[!muts_full$causal_sal=="causal"]))$auc.integral
  LEA3.2_lfmm2_AUCPR_sal_neutSNPs<- (pr.curve(scores.class0 =muts_full$LEA3.2_lfmm2_mlog10P_salenv[muts_full$causal_sal=="causal"], scores.class1 = muts_full$LEA3.2_lfmm2_mlog10P_salenv[muts_full$causal_sal=="neutral"]))$auc.integral
  # PRROC: computing and visualizing precision-recall and receiver operating characteristic curves in R
  # class 0 is the true positives, class 1 is the true negatives
  # it's important that the signal of class 0 is larger than class 1
  
  
## Correlate pop structure and allele frequencies ####
  # this takes a tick
  muts_full$structure_cor_G_PC1 <- muts_full$structure_cor_G_LFMM_U1_modtemp<- muts_full$structure_cor_G_LFMM_U1_modsal<- NA
  for (i in 1:nrow(muts_full)){
    muts_full$structure_cor_G_PC1[i] <- cor(G_full_subset[i,], subset_indPhen_df$PC1, method = "kendall")
    muts_full$structure_cor_G_LFMM_U1_modtemp[i] <- cor(G_full_subset[i,],mod_temp@U[,1], method = "kendall")
    muts_full$structure_cor_G_LFMM_U1_modsal[i] <- cor(G_full_subset[i,],mod_sal@U[,1], method = "kendall")  
  }
  
  if(!identical(as.character(muts_full$mutname), as.character(rownames(G_full_subset)))){print("Error 1i: mutations not lined up");break()}
  
  #VIZ:
  #sig by different methods vs. corr locus with structure
  #histogram of causal vs. neutral loci correlation with structure
  
  # cor plot of the different structure correlations

    subset_indPhen_df$LFMM_U1_temp <- mod_temp@U[,1]
    subset_indPhen_df$LFMM_U1_sal <- mod_sal@U[,1]
    
    if(K>1){
      subset_indPhen_df$LFMM_U2_temp <- mod_temp@U[,2]
      subset_indPhen_df$LFMM_U2_sal <- mod_sal@U[,2]
    }else{
      subset_indPhen_df$LFMM_U2_temp <- NA
      subset_indPhen_df$LFMM_U2_sal <- NA
    }
 

    cor_PC1_temp <- cor(subset_indPhen_df$PC1, subset_indPhen_df$temp_opt)
    cor_PC1_sal <- cor(subset_indPhen_df$PC1, subset_indPhen_df$sal_opt)
    cor_PC2_temp <- cor(subset_indPhen_df$PC2, subset_indPhen_df$temp_opt)
    cor_PC2_sal <- cor(subset_indPhen_df$PC2, subset_indPhen_df$sal_opt)
    cor_LFMMU1_temp <- cor(subset_indPhen_df$LFMM_U1_temp, subset_indPhen_df$temp_opt)
    cor_LFMMU1_sal <- cor(subset_indPhen_df$LFMM_U1_sal, subset_indPhen_df$sal_opt)
    cor_LFMMU2_temp <- cor(subset_indPhen_df$LFMM_U2_temp, subset_indPhen_df$temp_opt)
    cor_LFMMU2_sal <- cor(subset_indPhen_df$LFMM_U2_sal, subset_indPhen_df$sal_opt)
    cor_PC1_LFMMU1_temp <- cor(subset_indPhen_df$LFMM_U1_temp, subset_indPhen_df$PC1) 
    cor_PC1_LFMMU1_sal <- cor(subset_indPhen_df$LFMM_U1_sal, subset_indPhen_df$PC1) 
    cor_PC2_LFMMU1_temp <- cor(subset_indPhen_df$LFMM_U1_temp, subset_indPhen_df$PC2) 
    cor_PC2_LFMMU1_sal <- cor(subset_indPhen_df$LFMM_U1_sal, subset_indPhen_df$PC2) 
    
   pdf(paste0(path,seed,"_pdf_3bStructureCorrs.pdf"), width=15, height=15)
    
    # PC vs. Env
    k <- ggplot(subset_indPhen_df) + geom_point(aes(x=PC1, y=temp_opt, size=sal_opt), color="grey20") + geom_point(aes(x=PC1, y=temp_opt, color=temp_opt), size=2.5) + ylab("Deme temperature") + ggtheme + scale_colour_gradient2(high=rgb(1,0.4,0.2), low="cornflowerblue", mid=rgb(0.8,0.8,0.7), name="Temp") + labs(size="Env2") + ggtitle(paste0(plotmain))
    l <- ggplot(subset_indPhen_df) + geom_point(aes(x=PC1, y=sal_opt, size=sal_opt), color="grey20") + geom_point(aes(x=PC1, y=sal_opt, color=temp_opt), size=2.5) + ylab("Deme Env2") + ggtheme + scale_colour_gradient2(high=rgb(1,0.4,0.2), low="cornflowerblue", mid=rgb(0.8,0.8,0.7), name="Temp") + labs(size="Env2") + ggtitle(paste0(plotmain))
    
    # PC vs LFMM latent factors
    m <- ggplot(subset_indPhen_df) + geom_point(aes(x=PC1, y=LFMM_U1_temp,size=sal_opt), color="grey20") + geom_point(aes(x=PC1, y=LFMM_U1_temp, color=temp_opt), size=2.5) + ylab("Latent factor 1 LFMM temp model") + ggtheme + scale_colour_gradient2(high=rgb(1,0.4,0.2), low="cornflowerblue", mid=rgb(0.8,0.8,0.7), name="Temp") + labs(size="Env2") + ggtitle(paste0(plotmain))
    n <- ggplot(subset_indPhen_df) + geom_point(aes(x=PC1, y=LFMM_U1_sal,size=sal_opt), color="grey20") + geom_point(aes(x=PC1, y=LFMM_U1_sal, color=temp_opt), size=2.5) + ylab("Latent factor 1 LFMM Env2 model") + ggtheme + scale_colour_gradient2(high=rgb(1,0.4,0.2), low="cornflowerblue", mid=rgb(0.8,0.8,0.7), name="Temp") + labs(size="Env2") + ggtitle(paste0(plotmain))
    o <- ggplot(subset_indPhen_df) + geom_point(aes(x=PC2, y=LFMM_U1_temp,size=sal_opt), color="grey20") + geom_point(aes(x=PC2, y=LFMM_U1_temp, color=temp_opt), size=2.5) + ylab("Latent factor 1 LFMM temp model") + ggtheme + scale_colour_gradient2(high=rgb(1,0.4,0.2), low="cornflowerblue", mid=rgb(0.8,0.8,0.7), name="Temp") + labs(size="Env2") + ggtitle(paste0(plotmain))
    p <- ggplot(subset_indPhen_df) + geom_point(aes(x=PC2, y=LFMM_U1_sal,size=sal_opt), color="grey20") + geom_point(aes(x=PC2, y=LFMM_U1_sal, color=temp_opt), size=2.5) + ylab("Latent factor 1 LFMM Env2 model") + ggtheme + scale_colour_gradient2(high=rgb(1,0.4,0.2), low="cornflowerblue", mid=rgb(0.8,0.8,0.7), name="Temp") + labs(size="Env2") + ggtitle(paste0(plotmain))
    grid.arrange(k, l, m, n, o, p, nrow=3)
    
    m1 <- ggplot(subset_indPhen_df) + geom_point(aes(x=LFMM_U1_temp, y=temp_opt, size=sal_opt), color="grey20") + geom_point(aes(x=LFMM_U1_temp, y=temp_opt, color=temp_opt), size=2.5) + ylab("Deme temperature") + ggtheme + scale_colour_gradient2(high=rgb(1,0.4,0.2), low="cornflowerblue", mid=rgb(0.8,0.8,0.7), name="Temp") + labs(size="Env2") + ggtitle(paste0(plotmain))
    n1 <-ggplot(subset_indPhen_df) + geom_point(aes(x=LFMM_U1_sal, y=sal_opt, size=sal_opt), color="grey20") + geom_point(aes(x=LFMM_U1_sal, y=sal_opt, color=temp_opt), size=2.5) + ylab("Deme Env2") + ggtheme + scale_colour_gradient2(high=rgb(1,0.4,0.2), low="cornflowerblue", mid=rgb(0.8,0.8,0.7), name="Temp") + labs(size="Env2") + ggtitle(paste0(plotmain))
    o1 <- ggplot(subset_indPhen_df) + geom_point(aes(x=LFMM_U2_temp, y=temp_opt, size=sal_opt), color="grey20") + geom_point(aes(x=LFMM_U2_temp, y=temp_opt, color=temp_opt), size=2.5) + ylab("Deme temperature") + ggtheme + scale_colour_gradient2(high=rgb(1,0.4,0.2), low="cornflowerblue", mid=rgb(0.8,0.8,0.7), name="Temp") + labs(size="Env2") + ggtitle(paste0(plotmain))
    p1 <-  ggplot(subset_indPhen_df) + geom_point(aes(x=LFMM_U2_sal, y=temp_opt, size=sal_opt), color="grey20") + geom_point(aes(x=LFMM_U2_sal, y=sal_opt, color=temp_opt), size=2.5) + ylab("Deme Env2") + ggtheme + scale_colour_gradient2(high=rgb(1,0.4,0.2), low="cornflowerblue", mid=rgb(0.8,0.8,0.7), name="Temp") + labs(size="Env2") + ggtitle(paste0(plotmain))
    grid.arrange(k,l,m1,n1,o1,p1)
  
  # Plot Deme temperature of individual vs. PC1 loading of individual
  # Plot cor(Genotype, structure PC1) vs. -log10 P cor(af, temp)
  # Plot cor(Genotype, structure PC1) vs. -log10 P lfmm temp model

    # PC1 plots
    # Temp corrected
    r <- ggplot() + 
      geom_point(data=muts_full, aes(x=abs(structure_cor_G_PC1), y = LEA3.2_lfmm2_mlog10P_tempenv), color=muts_full$colors) +
      ggtheme +
      geom_point(data = muts_full[muts_full$LEA3.2_lfmm2_mlog10P_tempenv_sig,], aes(x=abs(structure_cor_G_PC1), y =LEA3.2_lfmm2_mlog10P_tempenv),pch=23, col=outlier_color, size=3, alpha=0.5) + 
      ylim(0, ymax)  + ylab("-log10(P) LFMM (Genotype, temp)") +
      ggtitle("Temperature - Structure corrected") + xlab("Abs(Cor(Genotype, Structure PC1))") + labs(color="temp VA prop.", size="temp VA prop.") +    
      geom_point(data = muts_full[muts_full$causal_temp=="causal",], aes(x=abs(structure_cor_G_PC1), y = LEA3.2_lfmm2_mlog10P_tempenv, 
                                                                         color=Va_temp_prop, size=Va_temp_prop), shape=shape_causal, alpha=0.5) + 
      scale_colour_viridis(option=cs, begin = begin_cs, end=end_cs, limits=c(0,1)) 
    
    # Temp uncorreected
    s <- ggplot() + 
      geom_point(data=muts_full, aes(x=abs(structure_cor_G_PC1), y = -log10(af_cor_temp_P)), color=muts_full$colors) + 
      ggtheme  +
      geom_point(data = muts_full[muts_full$cor_temp_sig,], aes(x=abs(structure_cor_G_PC1), y = -log10(af_cor_temp_P)),pch=23, col=adjustcolor("darkorange",0.5), size=3, alpha=0.5) + 
      ylim(0, ymax)  + ylab("-log10(P) cor(Genotype, temp)") +
      ggtitle("Temperature - Structure not corrected") + xlab("Abs(Cor(Genotype, Structure PC1))") + labs(color="temp VA prop.", size="temp VA prop.") +
      geom_point(data = muts_full[muts_full$causal_temp=="causal",], aes(x=abs(structure_cor_G_PC1), y = -log10(af_cor_temp_P), 
                                                                         color=Va_temp_prop, size=Va_temp_prop), shape=shape_causal, alpha=0.5) + 
      scale_colour_viridis(option=cs, begin = begin_cs, end=end_cs, limits=c(0,1))
   
    # Sal corrected
    rsal <- ggplot() + 
      geom_point(data=muts_full, aes(x=abs(structure_cor_G_PC1), y = LEA3.2_lfmm2_mlog10P_salenv), color=muts_full$colors) +
      ggtheme +
      geom_point(data = muts_full[muts_full$LEA3.2_lfmm2_mlog10P_salenv_sig,], aes(x=abs(structure_cor_G_PC1), y =LEA3.2_lfmm2_mlog10P_salenv),pch=23, col=outlier_color, size=3, alpha=0.5) + 
      ylim(0, ymax)  + ylab("-log10(P) LFMM (Genotype, Env2)") +
      ggtitle("Env2 - Structure corrected") + xlab("Abs(Cor(Genotype, Structure PC1))") + labs(color="Env2 VA prop.", size="Env2 VA prop.") +
      geom_point(data = muts_full[muts_full$causal_sal=="causal",], aes(x=abs(structure_cor_G_PC1), y = LEA3.2_lfmm2_mlog10P_salenv, 
                                                                        color=Va_sal_prop, size=Va_sal_prop), shape=shape_causal, alpha=0.5) + 
      scale_colour_viridis(option=cs, begin = begin_cs, end=end_cs, limits=c(0,1)) 
    
    # Sal uncorrected
    ssal <- ggplot() + 
      geom_point(data=muts_full, aes(x=abs(structure_cor_G_PC1), y = -log10(af_cor_sal_P)), color=muts_full$colors) + 
   ggtheme  +
       geom_point(data = muts_full[muts_full$cor_sal_sig,], aes(x=abs(structure_cor_G_PC1), y =-log10(af_cor_sal_P)),pch=23, col=adjustcolor("darkorange",0.5), size=3, alpha=0.5) + 
      ylim(0, ymax)  + ylab("-log10(P) cor(Genotype, Env2)") +
      ggtitle("Env2 - Structure not corrected") + xlab("Abs(Cor(Genotype, Structure PC1))") + labs(color="Env2 VA prop.", size="Env2 VA prop.") +
      geom_point(data = muts_full[muts_full$causal_sal=="causal",], aes(x=abs(structure_cor_G_PC1), y = -log10(af_cor_sal_P), 
                                                                        color=Va_sal_prop, size=Va_sal_prop), shape=shape_causal, alpha=0.5) + 
      scale_colour_viridis(option=cs, begin = begin_cs, end=end_cs, limits=c(0,1)) 
      
    
    
    grid.arrange(k, l, s, ssal, r, rsal, nrow=3)
    
    dev.off()
    
    print("LFMM finished")
  
    #cor.test(abs(muts_full$structure_cor_G_PC1), muts_full$LEA3.2_lfmm2_mlog10P_tempenv, method="pearson")
    
    # q <- ggplot() + 
    #   geom_point(data=muts_full, aes(x=abs(structure_cor_G_LFMM_U1_modtemp), y = LEA3.2_lfmm2_mlog10P_tempenv), color=muts_full$colors) +
    #   geom_point(data = muts_full[muts_full$causal_temp=="causal",], aes(x=abs(structure_cor_G_LFMM_U1_modtemp), y = LEA3.2_lfmm2_mlog10P_tempenv, 
    #                                                                      color=Va_temp_prop, size=Va_temp_prop), shape=shape_causal) + 
    #   scale_colour_viridis(option=cs, begin = begin_cs, end=end_cs, limits=c(0,1)) + 
    #   ggtheme +
    #   geom_point(data = muts_full[muts_full$LEA3.2_lfmm2_mlog10P_tempenv_sig,], aes(x=abs(structure_cor_G_LFMM_U1_modtemp), y =LEA3.2_lfmm2_mlog10P_tempenv),pch=23, col=outlier_color, size=3) + 
    #   ylim(0, ymax)  + ylab("-log10(P) LFMM temp") +
    #   ggtitle(paste0(plotmain," temp")) + xlab("Abs(Cor(Genotype, Structure LFMM U1 temp model))") + labs(color="temp VA prop.", size="temp VA prop.")
    # 
    # qaf <- ggplot() + 
    #   geom_point(data=muts_full, aes(x=abs(structure_cor_G_LFMM_U1_modtemp), y = -log10(af_cor_temp_P)), color=muts_full$colors) +
    #   geom_point(data = muts_full[muts_full$causal_temp=="causal",], aes(x=abs(structure_cor_G_LFMM_U1_modtemp), y = -log10(af_cor_temp_P), 
    #                                                                      color=Va_temp_prop, size=Va_temp_prop), shape=shape_causal) + 
    #   scale_colour_viridis(option=cs, begin = begin_cs, end=end_cs, limits=c(0,1)) + 
    #   ggtheme +
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
    #   ggtheme +
    #   ylim(0, ymax)  + ylab("-log10(P) LFMM temp") +
    #   ggtitle(paste0(plotmain," temp")) + xlab("-log10 P cor(af, temp)") + labs(color="temp VA prop.", size="temp VA prop.")
    # 
  
  ### LFMM Manhattan ####
  
  pdf(paste0(path,seed,"_pdf_4manhattan_LFMM.pdf"), width=15, height=10)
  
  cs <- "mako"
  begin_cs <- 0.9
  end_cs <- 0.2
  shape_causal <- 17
  outlier_color <- adjustcolor("brown1", 0.8)
  ymax=max(c(muts_full$LEA3.2_lfmm2_mlog10P_tempenv, muts_full$LEA3.2_lfmm2_mlog10P_salenv, 10))
  
  a <- ggplot() + 
    geom_point(data=muts_full, aes(x=pos_pyslim, y = LEA3.2_lfmm2_mlog10P_tempenv), color=muts_full$colors) + 
    ggtheme +
    geom_point(data = muts_full[muts_full$LEA3.2_lfmm2_mlog10P_tempenv_sig,], aes(x=pos_pyslim, y =LEA3.2_lfmm2_mlog10P_tempenv),pch=23, col=outlier_color, size=3, alpha=0.5) + 
    ylim(0, ymax)  + ylab("-log10(P) LFMM temp") +
    ggtitle(paste0(plotmain," temp")) + xlab("position") + labs(color="temp VA prop.", size="temp VA prop.") +
    geom_point(data = muts_full[muts_full$causal_temp=="causal",], aes(x=pos_pyslim, y = LEA3.2_lfmm2_mlog10P_tempenv, 
                                                                       color=Va_temp_prop, size=Va_temp_prop), shape=shape_causal, alpha=0.5) + 
    scale_colour_viridis(option=cs, begin = begin_cs, end=end_cs, limits=c(0,1))
  
  b <- ggplot() + 
    geom_point(data=muts_full, aes(x=pos_pyslim, y = LEA3.2_lfmm2_mlog10P_salenv), color=muts_full$colors)  + 
    ggtheme +
    geom_point(data = muts_full[muts_full$LEA3.2_lfmm2_mlog10P_salenv_sig,], aes(x=pos_pyslim, y =LEA3.2_lfmm2_mlog10P_salenv),pch=23, col=outlier_color, size=3, alpha=0.5) + 
    ylim(0, ymax)  + ylab("-log10(P) LFMM salinity") + ggtitle(paste0(plotmain," Env2")) + xlab("position") + labs(color="Env2 VA prop.", size="Env2 VA prop.") +
    geom_point(data = muts_full[muts_full$causal_sal=="causal",], aes(x=pos_pyslim, y = LEA3.2_lfmm2_mlog10P_salenv, 
                                                                      color=Va_sal_prop, size=Va_sal_prop), shape=shape_causal, alpha=0.5) + 
    scale_colour_viridis(option=cs, begin = begin_cs, end=end_cs, limits=c(0,1))
  
  grid.arrange(a, b, nrow=2)
  
  dev.off()
  
### RDA vegan no structure correction####
  
  sal <- subset_indPhen_df$sal_opt
  temp <- subset_indPhen_df$temp_opt
  
  rdaout <- rda(t(G_full_subset) ~ sal + temp)
  
  saveRDS(rdaout, paste0(path,seed,"_RDA.RDS"))
  
  
  #rdaout_corr <- rda(t(G_full_subset)~ sal + temp + Condition(subset_indPhen_df$x + subset_indPhen_df$y + subset_indPhen_df$PC1 + subset_indPhen_df$PC2)) #At first I tried adding geography (x,y location) and structure (PC axis) to the RDA, this completely corrected for everything, which led to low performance
  rdaout_corr <- rda(t(G_full_subset)~ sal + temp + Condition(subset_indPhen_df$PC1 + subset_indPhen_df$PC2)) #only structure correction
  
  saveRDS(rdaout_corr, paste0(path,seed,"_RDA_structcorr.RDS"))
  
  
  #str(rdaout)
  scores <- scores(rdaout, choices=1:4)
  scores_corr <- scores(rdaout_corr, choices=1:4)
  
  loci.sc <- scores$species
  loci.sc_corr <- scores_corr$species
  ind.sc <- scores$sites
  ind.sc_corr <- scores_corr$sites
  
  
  subset_indPhen_df$RDA1 <- ind.sc[,1]
  subset_indPhen_df$RDA2 <- ind.sc[,2]
  subset_indPhen_df$RDA1_corr <- ind.sc_corr[,1]
  subset_indPhen_df$RDA2_corr <- ind.sc_corr[,2]
  subset_indPhen_df$RDA_PC1 <-  ind.sc[,3]
  subset_indPhen_df$RDA_PC2 <-  ind.sc[,4]
  subset_indPhen_df$RDA_PC1_corr <-  ind.sc_corr[,3]
  subset_indPhen_df$RDA_PC2_corr <-  ind.sc_corr[,4]
  
  muts_full$RDA1_score <- loci.sc[,1]
  muts_full$RDA2_score <- loci.sc[,2]
  muts_full$RDA1_score_corr <- loci.sc_corr[,1]
  muts_full$RDA2_score_corr <- loci.sc_corr[,2]
  
  if(!identical(as.character(muts_full$mutname), as.character(rownames(G_full_subset)))){print("Error 1i: mutations not lined up");break()}
  
  
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
  ps_corr <- rdadapt(rdaout_corr, 2)
  
  muts_full$RDA_mlog10P <- -log10(ps$p.values)
  muts_full$RDA_mlog10P_sig <- ps$q.values<0.05
  muts_full$RDA_mlog10P_corr <- -log10(ps_corr$p.values)
  muts_full$RDA_mlog10P_sig_corr <- ps_corr$q.values<0.05
  
  # Compare p-vales for corrected vs. uncorrected RDA
  plot(muts_full$RDA_mlog10P, muts_full$RDA_mlog10P_corr)
  abline(0,1, col="red") #totally opposite outliers
  
  RDA_mlog10P_sig_noutliers <- sum(muts_full$RDA_mlog10P_sig)
  
  # Proportion VA
  (RDA_Va_temp_prop <- sum(muts_full$Va_temp_prop[muts_full$RDA_mlog10P_sig], na.rm=TRUE))
  (RDA_Va_temp_prop_corr <- sum(muts_full$Va_temp_prop[muts_full$RDA_mlog10P_sig_corr], na.rm=TRUE))
  
  (RDA_Va_sal_prop <- sum(muts_full$Va_sal_prop[muts_full$RDA_mlog10P_sig], na.rm=TRUE))
  (RDA_Va_sal_prop_corr <- sum(muts_full$Va_sal_prop[muts_full$RDA_mlog10P_sig_corr], na.rm=TRUE))
  
  (RDA_TPR <- sum(muts_full$RDA_mlog10P_sig & muts_full$causal)/sum(muts_full$causal))
  (RDA_TPR_corr <- sum(muts_full$RDA_mlog10P_sig_corr & muts_full$causal)/sum(muts_full$causal))
  
  (RDA_FDR_allSNPs <- sum(muts_full$RDA_mlog10P_sig & !muts_full$causal)/sum(muts_full$RDA_mlog10P_sig)) #non-causal outliers / (all outliers)
  (RDA_FDR_allSNPs_corr <- sum(muts_full$RDA_mlog10P_sig_cor & !muts_full$causal)/sum(muts_full$RDA_mlog10P_sig_corr)) #non-causal outliers / (all outliers)
  
  (num_RDA_sig_causal <- sum(muts_full$RDA_mlog10P_sig & muts_full$causal))
  (num_RDA_sig_neutral <- sum(muts_full$RDA_mlog10P_sig & muts_full$causal_temp=="neutral"))
 
  (num_RDA_sig_causal_corr <- sum(muts_full$RDA_mlog10P_sig_corr & muts_full$causal))
  (num_RDA_sig_neutral_corr <- sum(muts_full$RDA_mlog10P_sig_corr & muts_full$causal_temp=="neutral"))
  
  (RDA_FDR_neutSNPs <- num_RDA_sig_neutral/(num_RDA_sig_neutral + num_RDA_sig_causal)) #neutral outliers / (all outliers))
  (RDA_FDR_neutSNPs_corr <- num_RDA_sig_neutral_corr/(num_RDA_sig_neutral_corr + num_RDA_sig_causal_corr)) #neutral outliers / (all outliers)
  
  (RDA_FPR_neutSNPs <- num_RDA_sig_neutral/(sum(muts_full$causal_temp=="neutral")))
  (RDA_FPR_neutSNPs_corr <- num_RDA_sig_neutral_corr/(sum(muts_full$causal_temp=="neutral")))
  
  
  RDA_AUCPR_allSNPs<- (pr.curve(scores.class0 = muts_full$RDA_mlog10P[muts_full$causal], scores.class1 = muts_full$RDA_mlog10P[!muts_full$causal]))$auc.integral
  # class 0 is the true positives, class 1 is the true negatives
  # it's important that the signal of class 0 is larger than class 1
  
  (RDA_AUCPR_neutSNPs<- (pr.curve(scores.class0= muts_full$RDA_mlog10P[muts_full$causal], scores.class1 = muts_full$RDA_mlog10P[muts_full$causal_temp=="neutral"]))$auc.integral)
    # only neutral SNPs on 2nd half of genome not affected by selection
  (RDA_AUCPR_neutSNPs_corr<- (pr.curve(scores.class0= muts_full$RDA_mlog10P_corr[muts_full$causal], scores.class1 = muts_full$RDA_mlog10P_corr[muts_full$causal_temp=="neutral"]))$auc.integral)
  
  
# RDA outlier and individual plots ####
 

  pdf(paste0(path,seed,"_pdf_5RDA.pdf"), width=8, height=7)

    
  a<- screeplot(rdaout)
  str(a)
  a$y # save this it's the eigenvalues
  RDA1_propvar <- round(a$y[1]/sum(a$y),3)
  RDA2_propvar <- round(a$y[2]/sum(a$y),3)
  
  a_corr <- screeplot(rdaout_corr)
  a_corr$y # save this it's the eigenvalues
  RDA1_propvar_corr <- round(a_corr$y[1]/sum(a_corr$y),3)
  RDA2_propvar_corr <- round(a_corr$y[2]/sum(a_corr$y),3)
  
  ## Mutation RDA
  
  # ggplot() + 
  #   geom_point(data=muts_full, aes(x=RDA1_score, y = RDA2_score), color=adjustcolor("grey80", 0.2)) +
  #   geom_point(data = muts_full[muts_full$causal_temp=="causal",], 
  #              aes(x=RDA1_score, y = RDA2_score, color=Va_temp_prop, size=Va_temp_prop), shape=shape_causal) +
  #   scale_colour_viridis(option=cs, begin = begin_cs, end=end_cs, limits=c(0,1)) +  
  #   ggtheme +  geom_point(data = muts_full[muts_full$RDA_mlog10P_sig,], aes(x=RDA1_score, y = RDA2_score),
  #                                 pch=23, col=outlier_color, size=3) +  xlab("RDA 1") + ylab("RDA 2") + 
  #   ggtitle(paste0(plotmain, " mutations - temp")) + 
  #   geom_segment(data=arrows, aes(x=x, y=y, xend=dx, yend=dy),arrow=arrow(length = unit(0.2,"cm"))) + geom_text(data=arrows, aes(x=dx, y=dy, label=name), hjust="right", vjust="bottom")+ labs(color="temp VA prop.", size="temp VA prop.")
  # 
  
  
  
### RDA MUTATION PLOTS ####  
  lims <- max(abs(c(muts_full$mutSalEffect,muts_full$mutTempEffect)), na.rm=TRUE) #limits for mut effect size
  muts_scale <- max(c(abs(muts_full$RDA1_score),abs(muts_full$RDA2_score)))*1.3 #limits for the x and y axes
  # for nicer plot
  
  sal_arrow_muts = rdaout$CCA$biplot[which(rownames(rdaout$CCA$biplot)=="sal"),]
  temp_arrow_muts = rdaout$CCA$biplot[which(rownames(rdaout$CCA$biplot)=="temp"),]
  
  mult <- muts_scale/max(c(abs(sal_arrow_muts), abs(temp_arrow_muts)))
  

  arrows <- data.frame(x=rep(0,2), y = rep(0, 2), #start at origin
                       dx = c(sal_arrow_muts[1]*mult, temp_arrow_muts[1]*mult),
                       dy = c(sal_arrow_muts[2]*mult, temp_arrow_muts[2]*mult), name = c("Env2", "Temp"))
  
  mult2 <- c(max(abs(muts_full$RDA1_score))*1.3,max(abs(muts_full$RDA2_score))*1.3)/
    c(max(abs(rdaout$CCA$biplot[,1])), max(abs(rdaout$CCA$biplot[,2])))
  
  arrows2 <- data.frame(x=rep(0,2), y = rep(0, 2), #start at origin
                       dx = c(sal_arrow_muts[1]*mult2[1], temp_arrow_muts[1]*mult2[1]),
                       dy = c(sal_arrow_muts[2]*mult2[2], temp_arrow_muts[2]*mult2[2]), name = c("Env2", "Temp"))
  
  
  # Mutation effect in RDA space for temp
  ggplot() +
    geom_point(data=muts_full, aes(x=RDA1_score, y = RDA2_score), color=adjustcolor("grey80", 0.2)) +
    geom_point(data = muts_full[muts_full$causal_temp=="causal",], 
               aes(x=RDA1_score, y = RDA2_score, color=mutTempEffect, size=Va_temp_prop), shape=shape_causal) +
    scale_colour_viridis(option="turbo", #begin = begin_cs, end=end_cs, 
                         limits=c(-lims, lims)
                         )+  
    ggtheme +  geom_point(data = muts_full[muts_full$RDA_mlog10P_sig,], aes(x=RDA1_score, y = RDA2_score),
                                  pch=23, col=adjustcolor("black",0.3), size=3) +  
    xlab(paste0("RDA 1 (",  RDA1_propvar*100, "%)")) + ylab(paste0("RDA 2 (",  RDA2_propvar*100, "%)")) + 
    ggtitle(paste0(plotmain, " mutations - temp")) + 
    geom_segment(data=arrows2, aes(x=x, y=y, xend=dx, yend=dy),arrow=arrow(length = unit(0.2,"cm"))) + 
    geom_text(data=arrows2, aes(x=dx, y=dy, label=name), hjust="right", vjust="bottom")+ labs(color="mut Temp effect", size="temp VA prop.")
  
  ggplot() +
    geom_point(data=muts_full, aes(x=RDA1_score, y = RDA2_score), color=adjustcolor("grey80", 0.2)) +
    geom_point(data = muts_full[muts_full$causal_temp=="causal",], 
               aes(x=RDA1_score, y = RDA2_score, color=mutTempEffect, size=Va_temp_prop), shape=shape_causal) +
    scale_colour_viridis(option="turbo", #begin = begin_cs, end=end_cs, 
                         limits=c(-lims, lims)
    )+  
    ggtheme +  geom_point(data = muts_full[muts_full$RDA_mlog10P_sig,], aes(x=RDA1_score, y = RDA2_score),
                          pch=23, col=adjustcolor("black",0.3), size=3) +  
    xlab(paste0("RDA 1 (",  RDA1_propvar*100, "%)")) + ylab(paste0("RDA 2 (",  RDA2_propvar*100, "%)")) + 
    ggtitle(paste0(plotmain, " mutations - temp")) + 
    geom_segment(data=arrows, aes(x=x, y=y, xend=dx, yend=dy),arrow=arrow(length = unit(0.2,"cm"))) + 
    geom_text(data=arrows, aes(x=dx, y=dy, label=name), hjust="right", vjust="bottom")+ labs(color="mut Temp effect", size="temp VA prop.")
  
  # Mutation effect in corrected RDA space for temp
  
  muts_scale_corr <- max(c(abs(muts_full$RDA1_score_corr),abs(muts_full$RDA2_score_corr)))*1.3 #limits for the x and y axes
  # for nicer plot
  
  sal_arrow_muts_corr = rdaout_corr$CCA$biplot[which(rownames(rdaout_corr$CCA$biplot)=="sal"),]
  temp_arrow_muts_corr = rdaout_corr$CCA$biplot[which(rownames(rdaout_corr$CCA$biplot)=="temp"),]
  
  mult_corr <- muts_scale_corr/max(c(abs(sal_arrow_muts_corr), abs(temp_arrow_muts_corr)))
  
  arrows_corr <- data.frame(x=rep(0,2), y = rep(0, 2), #start at origin
                       dx = c(sal_arrow_muts_corr[1]*mult_corr, temp_arrow_muts_corr[1]*mult_corr),
                       dy = c(sal_arrow_muts_corr[2]*mult_corr, temp_arrow_muts_corr[2]*mult_corr), name = c("Env2", "Temp"))
  
  mult2_corr <- c(max(abs(muts_full$RDA1_score_corr))*1.3,max(abs(muts_full$RDA2_score_corr))*1.3)/
    c(max(abs(rdaout_corr$CCA$biplot[,1])), max(abs(rdaout_corr$CCA$biplot[,2])))
  
  arrows2_corr <- data.frame(x=rep(0,2), y = rep(0, 2), #start at origin
                        dx = c(sal_arrow_muts_corr[1]*mult2_corr[1], temp_arrow_muts_corr[1]*mult2_corr[1]),
                        dy = c(sal_arrow_muts_corr[2]*mult2_corr[2], temp_arrow_muts_corr[2]*mult2_corr[2]), name = c("Env2", "Temp"))
  
  
  
  ggplot() +
    geom_point(data=muts_full, aes(x=RDA1_score_corr, y = RDA2_score_corr), color=adjustcolor("grey80", 0.2)) +
    geom_point(data = muts_full[muts_full$causal_temp=="causal",], 
               aes(x=RDA1_score_corr, y = RDA2_score_corr, color=mutTempEffect, size=Va_temp_prop), shape=shape_causal) +
    scale_colour_viridis(option="turbo", #begin = begin_cs, end=end_cs, 
                         limits=c(-lims, lims)
    )+  
    ggtheme +  geom_point(data = muts_full[muts_full$RDA_mlog10P_sig_corr,], aes(x=RDA1_score_corr, y = RDA2_score_corr),
                          pch=23, col=adjustcolor("black",0.3), size=3) +  
    xlab(paste0("RDA 1 (",  RDA1_propvar_corr*100, "%)")) + ylab(paste0("RDA 2 (",  RDA2_propvar_corr*100, "%)")) + 
    ggtitle(paste0(plotmain, " mutations - temp\nCorrected RDA")) + 
    geom_segment(data=arrows_corr, aes(x=x, y=y, xend=dx, yend=dy),arrow=arrow(length = unit(0.2,"cm"))) + 
    geom_text(data=arrows_corr, aes(x=dx, y=dy, label=name), hjust="right", vjust="bottom")+ labs(color="mut Temp effect", size="temp VA prop.")
  
  ggplot() +
    geom_point(data=muts_full, aes(x=RDA1_score_corr, y = RDA2_score_corr), color=adjustcolor("grey80", 0.2)) +
    geom_point(data = muts_full[muts_full$causal_temp=="causal",], 
               aes(x=RDA1_score_corr, y = RDA2_score_corr, color=mutTempEffect, size=Va_temp_prop), shape=shape_causal) +
    scale_colour_viridis(option="turbo", #begin = begin_cs, end=end_cs, 
                         limits=c(-lims, lims)
    )+  
    ggtheme +  geom_point(data = muts_full[muts_full$RDA_mlog10P_sig_corr,], aes(x=RDA1_score_corr, y = RDA2_score_corr),
                          pch=23, col=adjustcolor("black",0.3), size=3) +  
    xlab(paste0("RDA 1 (",  RDA1_propvar_corr*100, "%)")) + ylab(paste0("RDA 2 (",  RDA2_propvar_corr*100, "%)")) + 
    ggtitle(paste0(plotmain, " mutations - temp\nCorrected RDA")) + 
    geom_segment(data=arrows2_corr, aes(x=x, y=y, xend=dx, yend=dy),arrow=arrow(length = unit(0.2,"cm"))) + 
    geom_text(data=arrows2_corr, aes(x=dx, y=dy, label=name), hjust="right", vjust="bottom")+ labs(color="mut Temp effect", size="temp VA prop.")
  
  
  ## Mutation effect in RDA space for salinity
  ggplot() + 
    geom_point(data=muts_full, aes(x=RDA1_score, y = RDA2_score), color=adjustcolor("grey80", 0.2)) +
    geom_point(data = muts_full[muts_full$causal_sal=="causal",], 
               aes(x=RDA1_score, y = RDA2_score, color=mutSalEffect, size=Va_sal_prop), shape=shape_causal) +
    scale_colour_viridis(option="turbo", #begin = begin_cs, end=end_cs, 
                         limits=c(-lims, lims)
    )+  
    ggtheme +  geom_point(data = muts_full[muts_full$RDA_mlog10P_sig,], aes(x=RDA1_score, y = RDA2_score),
                                  pch=23, col=adjustcolor("black",0.3), size=3) +  
    xlab(paste0("RDA 1 (",  RDA1_propvar*100, "%)")) + ylab(paste0("RDA 2 (",  RDA2_propvar*100, "%)")) + 
    ggtitle(paste0(plotmain, " mutations - Env2")) + 
    geom_segment(data=arrows2, aes(x=x, y=y, xend=dx, yend=dy),arrow=arrow(length = unit(0.2,"cm"))) + 
    geom_text(data=arrows2, aes(x=dx, y=dy, label=name), hjust="right", vjust="bottom")+ 
    labs(color="mut Env2 effect", size="Env2 VA prop.")
  
  ## Mutation effect in corrected RDA space for salinity
  ggplot() + 
    geom_point(data=muts_full, aes(x=RDA1_score_corr, y = RDA2_score_corr), color=adjustcolor("grey80", 0.2)) +
    geom_point(data = muts_full[muts_full$causal_sal=="causal",], 
               aes(x=RDA1_score_corr, y = RDA2_score_corr, color=mutSalEffect, size=Va_sal_prop), shape=shape_causal) +
    scale_colour_viridis(option="turbo", #begin = begin_cs, end=end_cs, 
                         limits=c(-lims, lims)
    )+  
    ggtheme +  geom_point(data = muts_full[muts_full$RDA_mlog10P_sig_corr,], aes(x=RDA1_score_corr, y = RDA2_score_corr),
                          pch=23, col=adjustcolor("black",0.3), size=3) +  
    xlab(paste0("RDA 1 (",  RDA1_propvar_corr*100, "%)")) + ylab(paste0("RDA 2 (",  RDA2_propvar_corr*100, "%)")) + 
    ggtitle(paste0(plotmain, " mutations - Env2\nCorrected RDA")) + 
    geom_segment(data=arrows2_corr, aes(x=x, y=y, xend=dx, yend=dy),arrow=arrow(length = unit(0.2,"cm"))) + 
    geom_text(data=arrows2_corr, aes(x=dx, y=dy, label=name), hjust="right", vjust="bottom")+ 
    labs(color="mut Env2 effect", size="Env2 VA prop.")
  
## RDA individual plot ####
  ind_scale <- max(c(abs(subset_indPhen_df$RDA1)), max(abs(subset_indPhen_df$RDA2)))
  # for nicer plot
  
  sal_arrow_muts = rdaout$CCA$biplot[which(rownames(rdaout$CCA$biplot)=="sal"),]
  temp_arrow_muts = rdaout$CCA$biplot[which(rownames(rdaout$CCA$biplot)=="temp"),]
  
  mult <- ind_scale/max(c(abs(sal_arrow_muts), abs(temp_arrow_muts)))
  
  arrows_ind <- data.frame(x=rep(0,2), y = rep(0, 2),
                           dx = c(sal_arrow_muts[1]*mult, temp_arrow_muts[1]*mult),
                           dy = c(sal_arrow_muts[2]*mult, temp_arrow_muts[2]*mult), name = c("Env2", "Temp"))
  
  mult2_ind<- c(max(abs(subset_indPhen_df$RDA1))*1.3,max(abs(subset_indPhen_df$RDA2))*1.3)/
    c(max(abs(rdaout$CCA$biplot[,1])), max(abs(rdaout$CCA$biplot[,2])))
  
  arrows2_ind <- data.frame(x=rep(0,2), y = rep(0, 2), #start at origin
                             dx = c(sal_arrow_muts[1]*mult2_ind[1], temp_arrow_muts[1]*mult2_ind[1]),
                             dy = c(sal_arrow_muts[2]*mult2_ind[2], temp_arrow_muts[2]*mult2_ind[2]), name = c("Env2", "Temp"))
  
  # Individuals in RDA space
  ggplot(subset_indPhen_df) + 
    geom_point(aes(x=RDA1, y = RDA2, size=sal_opt), color="grey20") + 
    geom_point(aes(x=RDA1, y=RDA2, color=temp_opt), size=2.5) + ggtheme + 
    xlab(paste0("RDA 1 (",  RDA1_propvar*100, "%)")) + ylab(paste0("RDA 2 (",  RDA2_propvar*100, "%)")) +
    scale_colour_gradient2(high=rgb(1,0.4,0.2), low="cornflowerblue", mid=rgb(0.8,0.8,0.7), name="Temp") + 
    labs(size="Env2") + 
    ggtitle(paste0(plotmain,"; Individs.; N traits = ", info$Ntraits)) + 
    geom_segment(data=arrows_ind, aes(x=x, y=y, xend=dx, yend=dy),arrow=arrow(length = unit(0.2,"cm"))) + 
    geom_text(data=arrows_ind, aes(x=dx, y=dy, label=name), hjust="right", vjust="bottom")
  
  ggplot(subset_indPhen_df) + 
    geom_point(aes(x=RDA1, y = RDA2, size=sal_opt), color="grey20") + 
    geom_point(aes(x=RDA1, y=RDA2, color=temp_opt), size=2.5) + ggtheme + 
    xlab(paste0("RDA 1 (",  RDA1_propvar*100, "%)")) + ylab(paste0("RDA 2 (",  RDA2_propvar*100, "%)")) +
    scale_colour_gradient2(high=rgb(1,0.4,0.2), low="cornflowerblue", mid=rgb(0.8,0.8,0.7), name="Temp") + 
    labs(size="Env2") + 
    ggtitle(paste0(plotmain,"; Individs.; N traits = ", info$Ntraits)) + 
    geom_segment(data=arrows2_ind, aes(x=x, y=y, xend=dx, yend=dy),arrow=arrow(length = unit(0.2,"cm"))) + 
    geom_text(data=arrows2_ind, aes(x=dx, y=dy, label=name), hjust="right", vjust="bottom")
  
  
  # Individuals in corrected RDA space
  
  ind_scale_corr <- max(c(abs(subset_indPhen_df$RDA1_corr), abs(subset_indPhen_df$RDA2_corr)))
  # for nicer plot
  sal_arrow_muts_corr = rdaout_corr$CCA$biplot[which(rownames(rdaout_corr$CCA$biplot)=="sal"),]
  temp_arrow_muts_corr = rdaout_corr$CCA$biplot[which(rownames(rdaout_corr$CCA$biplot)=="temp"),]
  mult_corr <- ind_scale_corr/max(c(abs(sal_arrow_muts_corr), abs(temp_arrow_muts_corr)))
  
  arrows_ind_corr <- data.frame(x=rep(0,2), y = rep(0, 2),
                                dx = c(sal_arrow_muts_corr[1]*mult_corr, temp_arrow_muts_corr[1]*mult_corr),
                                dy = c(sal_arrow_muts_corr[2]*mult_corr, temp_arrow_muts_corr[2]*mult_corr), name = c("Env2", "Temp"))
  
  mult2_ind_corr<- c(max(abs(subset_indPhen_df$RDA1_corr))*1.3,max(abs(subset_indPhen_df$RDA2_corr))*1.3)/
    c(max(abs(rdaout_corr$CCA$biplot[,1])), max(abs(rdaout_corr$CCA$biplot[,2])))
  
  arrows2_ind_corr <- data.frame(x=rep(0,2), y = rep(0, 2), #start at origin
                            dx = c(sal_arrow_muts_corr[1]*mult2_ind_corr[1], temp_arrow_muts_corr[1]*mult2_ind_corr[1]),
                            dy = c(sal_arrow_muts_corr[2]*mult2_ind_corr[2], temp_arrow_muts_corr[2]*mult2_ind_corr[2]), name = c("Env2", "Temp"))
  
  
  ggplot(subset_indPhen_df) + 
    geom_point(aes(x=RDA1_corr, y = RDA2_corr, size=sal_opt), color="grey20") + 
    geom_point(aes(x=RDA1_corr, y=RDA2_corr, color=temp_opt), size=2.5) + ggtheme + 
    xlab(paste0("RDA 1 (",  RDA1_propvar_corr*100, "%)")) + ylab(paste0("RDA 2 (",  RDA2_propvar_corr*100, "%)")) +
    scale_colour_gradient2(high=rgb(1,0.4,0.2), low="cornflowerblue", mid=rgb(0.8,0.8,0.7), name="Temp") + 
    labs(size="Env2") + 
    ggtitle(paste0(plotmain,"; Individs.; N traits = ", info$Ntraits,"\nRDA corrected")) + 
    geom_segment(data=arrows_ind_corr, aes(x=x, y=y, xend=dx, yend=dy),arrow=arrow(length = unit(0.2,"cm"))) + 
    geom_text(data=arrows_ind_corr, aes(x=dx, y=dy, label=name), hjust="right", vjust="bottom")
  
  ggplot(subset_indPhen_df) + 
    geom_point(aes(x=RDA1_corr, y = RDA2_corr, size=sal_opt), color="grey20") + 
    geom_point(aes(x=RDA1_corr, y=RDA2_corr, color=temp_opt), size=2.5) + ggtheme + 
    xlab(paste0("RDA 1 (",  RDA1_propvar_corr*100, "%)")) + ylab(paste0("RDA 2 (",  RDA2_propvar_corr*100, "%)")) +
    scale_colour_gradient2(high=rgb(1,0.4,0.2), low="cornflowerblue", mid=rgb(0.8,0.8,0.7), name="Temp") + 
    labs(size="Env2") + 
    ggtitle(paste0(plotmain,"; Individs.; N traits = ", info$Ntraits,"\nRDA corrected")) + 
    geom_segment(data=arrows2_ind_corr, aes(x=x, y=y, xend=dx, yend=dy),arrow=arrow(length = unit(0.2,"cm"))) + 
    geom_text(data=arrows2_ind_corr, aes(x=dx, y=dy, label=name), hjust="right", vjust="bottom")
  
  dev.off()
  
#### RDA manhattan ####
  ymax = max(c(muts_full$RDA_mlog10P, 10))
  pdf(paste0(path,seed,"_pdf_6manhattan_RDA.pdf"), width=15, height=15)
  
  a<- ggplot() + 
    geom_point(data=muts_full, 
               aes(x=pos_pyslim, y = RDA_mlog10P), color=muts_full$colors) +
    geom_point(data = muts_full[muts_full$RDA_mlog10P_sig,], aes(x=pos_pyslim, y =RDA_mlog10P),pch=23, col="grey", size=3, alpha=0.5) + 
    geom_point(data = muts_full[muts_full$causal_temp=="causal",], 
               aes(x=pos_pyslim, y = RDA_mlog10P, color= mutTempEffect, 
                   size=Va_temp_prop), shape=shape_causal) + 
    scale_colour_viridis(option="turbo") + 
    ggtheme +
    ylim(0, ymax)  + ylab("-log10(P) RDA") + ggtitle(paste0(plotmain," temp")) + xlab("position") + labs(color="MutTempEffect", size="Temp VA prop.")
  
  b<- ggplot() + 
    geom_point(data=muts_full, aes(x=pos_pyslim, y = RDA_mlog10P), color=muts_full$colors) +
    geom_point(data = muts_full[muts_full$RDA_mlog10P_sig,], aes(x=pos_pyslim, y =RDA_mlog10P),pch=23, col="grey", size=3, alpha=0.5) + 
    geom_point(data = muts_full[muts_full$causal_sal=="causal",], 
               aes(x=pos_pyslim, y = RDA_mlog10P, color=mutSalEffect, size=Va_sal_prop), 
               shape=shape_causal) + 
    scale_colour_viridis(option="turbo") + 
    ggtheme +
    ylim(0, ymax)  + ylab("-log10(P) RDA") + ggtitle(paste0(plotmain," Env2")) + xlab("position") + labs(color="MutSalEffect", size="Env2 VA prop.")
  
  c<- ggplot() + 
    geom_point(data=muts_full, 
               aes(x=pos_pyslim, y = RDA_mlog10P_corr), color=muts_full$colors) +
    geom_point(data = muts_full[muts_full$RDA_mlog10P_sig_corr,], aes(x=pos_pyslim, y =RDA_mlog10P_corr),pch=23, col="grey", size=3, alpha=0.5) + 
    geom_point(data = muts_full[muts_full$causal_temp=="causal",], 
               aes(x=pos_pyslim, y = RDA_mlog10P_corr, color= mutTempEffect, 
                   size=Va_temp_prop), shape=shape_causal) + 
    scale_colour_viridis(option="turbo") + 
    ggtheme +
    ylim(0, ymax)  + ylab("-log10(P) RDA") + ggtitle(paste0(plotmain," temp\n RDA corrected")) + xlab("position") + labs(color="MutTempEffect", size="Temp VA prop.")
  
  d<- ggplot() + 
    geom_point(data=muts_full, aes(x=pos_pyslim, y = RDA_mlog10P_corr), color=muts_full$colors) +
    geom_point(data = muts_full[muts_full$RDA_mlog10P_sig_corr,], aes(x=pos_pyslim, y =RDA_mlog10P_corr),pch=23, col="grey", size=3, alpha=0.5) + 
    geom_point(data = muts_full[muts_full$causal_sal=="causal",], 
               aes(x=pos_pyslim, y = RDA_mlog10P_corr, color=mutSalEffect, size=Va_sal_prop), 
               shape=shape_causal) + 
    scale_colour_viridis(option="turbo") + 
    ggtheme +
    ylim(0, ymax)  + ylab("-log10(P) RDA") + ggtitle(paste0(plotmain," Env2\nRDA corrected")) + xlab("position") + labs(color="MutSalEffect", size="Env2 VA prop.")
  
  grid.arrange(a, b,c,d, nrow=4)
  dev.off()
  
#### write RDA table ####
  write.table(data.frame(seed=as.character(seed), 
                         summary(rdaout)$biplot),paste0(path,seed,"_RDAloadings.txt")) 
  
  write.table(data.frame(seed=as.character(seed), 
                         summary(rdaout_corr)$biplot),paste0(path,seed,"_RDAcorrected_loadings.txt")) 
  
  
  if(!identical(as.character(muts_full$mutname), as.character(rownames(G_full_subset)))){print("Error 1j: mutations not lined up");break()}
  
### RDA predict phenotype from random loci ####
  num_neut = num_neut_postfilter
  nmax <- min(c(20000, num_causal_postfilter+num_neut))
  
  if ((num_causal_postfilter+num_neut)>4999){
    nloci <- c(10, 50, 100, 500, 5000, nmax)
  }else{
    nloci <- c(10, 50, 100, 500, nmax)
  }
  
  
  Random_cor_RDAsalpredict_salphen <- Random_cor_RDAtemppredict_tempphen_structcorr <- rep(NA, length(nloci))
  Random_cor_RDAtemppredict_tempphen <- Random_cor_RDAsalpredict_salphen_structcorr <- rep(NA, length(nloci))
  
  for (i in 1:length(nloci)){
    
    ## No structure correction
    whichlocs <- sort(sample(1:nrow(muts_full), nloci[i], replace=FALSE))
    G_random_out <- G_full_subset[sort(whichlocs),]
    rdaout_rand <- rda(t(G_random_out)~ subset_indPhen_df$sal_opt + subset_indPhen_df$temp_opt)
    scores_rand <- scores(rdaout_rand, choices=1:4)
    loci.sc_rand <- scores_rand$species
    ind.sc_rand <- scores_rand$sites

    # if (nloci[i]==5000){ # make a plot
    # subset_indPhen_df$RDA1_rand <- ind.sc_rand[,1]
    # subset_indPhen_df$RDA2_rand <- ind.sc_rand[,2]
    #   a_rand <- screeplot(rdaout_rand)
    #   a_rand$y # save this it's the eigenvalues
    #   RDA1_propvar_rand <- round(a_rand$y[1]/sum(a_rand$y),3)
    #   RDA2_propvar_rand <- round(a_rand$y[2]/sum(a_rand$y),3)
    #   sal_arrow_muts = rdaout_rand$CCA$biplot[grep("sal",rownames(rdaout_rand$CCA$biplot)),]
    #   temp_arrow_muts = rdaout_rand$CCA$biplot[grep("temp",rownames(rdaout_rand$CCA$biplot)),]
    #   arrows_ind <- data.frame(x=rep(0,2), y = rep(0, 2),
    #                            dx = c(sal_arrow_muts[1], temp_arrow_muts[1]),
    #                            dy = c(sal_arrow_muts[2], temp_arrow_muts[2]), name = c("Env2", "Temp"))
    # 
    #   ggplot(subset_indPhen_df) + 
    #     geom_point(aes(x=RDA1_rand, y = RDA2_rand, size=sal_opt), color="grey20") + 
    #     geom_point(aes(x=RDA1_rand, y=RDA2_rand, color=temp_opt), size=2.5) + ggtheme + 
    #     xlab(paste0("RDA 1 (",  RDA1_propvar_rand*100, "%)")) + ylab(paste0("RDA 2 (",  RDA2_propvar_rand*100, "%)")) +
    #     scale_colour_gradient2(high=rgb(1,0.4,0.2), low="cornflowerblue", mid=rgb(0.8,0.8,0.7), name="Temp") + 
    #     labs(size="Env2") + 
    #     ggtitle(paste0(plotmain,"; Individs.; N traits = ", info$Ntraits, "\nRDA random loci: ",nloci[i])) + 
    #     geom_segment(data=arrows_ind, aes(x=x, y=y, xend=dx, yend=dy),arrow=arrow(length = unit(0.2,"cm"))) + 
    #     geom_text(data=arrows_ind, aes(x=dx, y=dy, label=name), hjust="right", vjust="bottom")
    #   
    # }
   
    
    temp_pred <- scale(ind.sc_rand[,1]*eigenvals(rdaout_rand)[1]*summary(rdaout_rand)$biplot[2,1] + 
                  ind.sc_rand[,2]*eigenvals(rdaout_rand)[2]*summary(rdaout_rand)$biplot[2,2])
    sal_pred <- scale(ind.sc_rand[,1]*eigenvals(rdaout_rand)[1]*summary(rdaout_rand)$biplot[1,1] + 
                ind.sc_rand[,2]*eigenvals(rdaout_rand)[2]*summary(rdaout_rand)$biplot[1,2])
    
    Random_cor_RDAtemppredict_tempphen[i] <- cor(scale(subset_indPhen_df$phen_temp), scale(temp_pred), method = "kendall") #scaling doesn't matter for this calculation
    Random_cor_RDAsalpredict_salphen[i] <- cor(scale(subset_indPhen_df$phen_sal), scale(sal_pred), method = "kendall")
    
    ## With structure correction
    rdaout_rand_corr <- rda(t(G_random_out)~ subset_indPhen_df$sal_opt + subset_indPhen_df$temp_opt + Condition(subset_indPhen_df$PC1 + subset_indPhen_df$PC2))
    scores_rand_corr <- scores(rdaout_rand_corr, choices=1:4)
    loci.sc_rand_corr <- scores_rand_corr$species
    ind.sc_rand_corr <- scores_rand_corr$sites
    
    # The following works because salinity is always 1st row of biplot and temp is always 2nd row - hard coding
    temp_pred_corr <- scale(ind.sc_rand_corr[,1]*eigenvals(rdaout_rand_corr)[1]*summary(rdaout_rand_corr)$biplot[2,1] + 
      ind.sc_rand_corr[,2]*eigenvals(rdaout_rand_corr)[2]*summary(rdaout_rand_corr)$biplot[2,2])
    sal_pred_corr <- scale(ind.sc_rand_corr[,1]*eigenvals(rdaout_rand_corr)[1]*summary(rdaout_rand_corr)$biplot[1,1] + 
      ind.sc_rand_corr[,2]*eigenvals(rdaout_rand_corr)[2]*summary(rdaout_rand_corr)$biplot[1,2])
    
    Random_cor_RDAtemppredict_tempphen_structcorr[i] <- cor(scale(subset_indPhen_df$phen_temp), scale(temp_pred_corr), method = "kendall")
    Random_cor_RDAsalpredict_salphen_structcorr[i] <- cor(scale(subset_indPhen_df$phen_sal), scale(sal_pred_corr), method = "kendall")
  }
  
  RDA1_temp_cor <- summary(rdaout)$biplot[2,1]
  RDA1_sal_cor <- summary(rdaout)$biplot[1,1]
  RDA2_temp_cor <- summary(rdaout)$biplot[2,2]
  RDA2_sal_cor <- summary(rdaout)$biplot[1,2]
  
  subset_indPhen_df$RDA_predict_tempPhen_20KSNPs <- scale(temp_pred)
  subset_indPhen_df$RDA_predict_salPhen_20KSNPs <- scale(sal_pred)
  subset_indPhen_df$RDA_predict_tempPhen_20KSNPs_structcorr <- scale(temp_pred_corr)
  subset_indPhen_df$RDA_predict_salPhen_20KSNPs_structcorr <- scale(sal_pred_corr)
  
  out_RDA_pred <- data.frame(nloci, Random_cor_RDAtemppredict_tempphen, 
                             Random_cor_RDAsalpredict_salphen,
                             Random_cor_RDAtemppredict_tempphen_structcorr, 
                             Random_cor_RDAsalpredict_salphen_structcorr)
  write.table(out_RDA_pred,paste0(path,seed,"_Rout_RDA_predictions.txt"))
  
  # Save the summary stat for the performace of the RDA prediction
  RDA_cor_RDA20000temppredict_tempPhen <- Random_cor_RDAtemppredict_tempphen[nloci==nmax]
  RDA_cor_RDA20000salpredict_salPhen <- Random_cor_RDAsalpredict_salphen[nloci==nmax]
  RDA_cor_RDA20000temppredict_tempPhen_structcorr <- Random_cor_RDAtemppredict_tempphen_structcorr[nloci==nmax]
  RDA_cor_RDA20000salpredict_salPhen_structcorr <- Random_cor_RDAsalpredict_salphen_structcorr[nloci==nmax]  
  
  ## RDA mutation prediction, correlation with mutation effect and Va ####
  
  # Without structure correction, using the RDA run on all the SNPs (rdaout object)
    # sum(locus_score * eigenvalue * loading_env_RDA_axis)
    # The following works because salinity is always 1st row of biplot and temp is always 2nd row - hard coding
    muts_full$RDA_mut_temp_pred <- scale(muts_full$RDA1_score*eigenvals(rdaout)[1]*summary(rdaout)$biplot[2,1] + #RDA1
                                   muts_full$RDA2_score*eigenvals(rdaout)[2]*summary(rdaout)$biplot[2,2])   #RDA2
    muts_full$RDA_mut_sal_pred <- scale(muts_full$RDA1_score*eigenvals(rdaout)[1]*summary(rdaout)$biplot[1,1] +
                                  muts_full$RDA2_score*eigenvals(rdaout)[2]*summary(rdaout)$biplot[1,2]) 
    #plot(muts_full$mutTempEffect, muts_full$RDA1_score)
    RDA_RDAmutpred_cor_tempEffect <- cor(muts_full$RDA_mut_temp_pred, muts_full$mutTempEffect, use="complete.obs", method="kendall")
    RDA_absRDAmutpred_cor_tempVa <- cor(abs(muts_full$RDA_mut_temp_pred), muts_full$Va_temp, use="complete.obs", method="kendall")
    
    if (info$Ntraits==1){
      RDA_RDAmutpred_cor_salEffect <- NA
      RDA_absRDAmutpred_cor_salVa <- NA
      }else{
      RDA_RDAmutpred_cor_salEffect <- cor(muts_full$RDA_mut_sal_pred, muts_full$mutSalEffect, use="complete.obs", method="kendall")
      RDA_absRDAmutpred_cor_salVa <- cor(abs(muts_full$RDA_mut_sal_pred), muts_full$Va_sal, use="complete.obs", method="kendall")
      }
    
  # With structure correction, using the RDA run on all the SNPs (rdaout_corr object)
    # sum(locus_score * eigenvalue * loading_env_RDA_axis)
    # The following works because salinity is always 1st row of biplot and temp is always 2nd row - hard coding
    muts_full$RDA_mut_temp_pred_structcorr <- scale(muts_full$RDA1_score_corr*eigenvals(rdaout_corr)[1]*summary(rdaout_corr)$biplot[2,1] + #RDA1
                                              muts_full$RDA2_score_corr*eigenvals(rdaout_corr)[2]*summary(rdaout_corr)$biplot[2,2])   #RDA2
    muts_full$RDA_mut_sal_pred_structcorr <- scale(muts_full$RDA1_score_corr*eigenvals(rdaout_corr)[1]*summary(rdaout_corr)$biplot[1,1] +
                                              muts_full$RDA2_score_corr*eigenvals(rdaout_corr)[2]*summary(rdaout_corr)$biplot[1,2]) 
    #plot(muts_full$mutTempEffect, muts_full$RDA1_score)
    
    RDA_RDAmutpred_cor_tempEffect_structcorr <- cor(muts_full$RDA_mut_temp_pred_structcorr, muts_full$mutTempEffect, use="complete.obs")
    RDA_absRDAmutpred_cor_tempVa_structcorr <- cor(abs(muts_full$RDA_mut_temp_pred_structcorr), muts_full$Va_temp, use="complete.obs")
    
    if (info$Ntraits==1){
      RDA_RDAmutpred_cor_salEffect_structcorr <- NA
      RDA_absRDAmutpred_cor_salVa_structcorr <- NA
    }else{
      RDA_RDAmutpred_cor_salEffect_structcorr <- cor(muts_full$RDA_mut_sal_pred_structcorr, muts_full$mutSalEffect, use="complete.obs")
      RDA_absRDAmutpred_cor_salVa_structcorr <- cor(abs(muts_full$RDA_mut_sal_pred_structcorr), muts_full$Va_sal, use="complete.obs")
    }
    
  pdf(paste0(path,seed,"_pdf_6zRDA_predict.pdf"), width=7, height=7)
  par(mfrow=c(1,1))
  
  plot(nloci, Random_cor_RDAtemppredict_tempphen, type="l", ylim=c(-1,1), ylab="Cor(RDA prediction, true phenotype)", bty="l" , 
       xlab="Number random loci in RDA", col="darkred", lwd=4, xlim=c(0, 20000), main = paste0(plotmain), cex.main=0.5, las=2)
  points(nloci, Random_cor_RDAsalpredict_salphen, type="l", col="cornflowerblue", lwd=2, lty=1)
  
  points(nloci, Random_cor_RDAtemppredict_tempphen_structcorr, type="l", col="darkred", lwd=2, lty=2)
  points(nloci, Random_cor_RDAsalpredict_salphen_structcorr, type="l", col="cornflowerblue", lwd=2, lty=2)
  
  legend(0, -0.25, c("Temp", "Env2"), col=c("darkred", "cornflowerblue"), lty = c(1,1), bty="l", lwd=c(4,1))
  legend(10000, -0.25, c("No structure correction", "Structure correction"), col="grey", lty = c(1,2), bty="l", lwd=2)
  
  plot(scale(muts_full$mutTempEffect), scale(muts_full$RDA_mut_temp_pred), main = paste0(plotmain), cex.main=0.5)
  abline(0,1, col="red")
  
  if (info$Ntraits>1){
    plot(scale(muts_full$mutSalEffect), scale(muts_full$RDA_mut_sal_pred), main =  paste0(plotmain), cex.main=0.5)
    abline(0,1, col="red")
  }
  
  dev.off()
  

### GWAS ####  
  dim(subset_indPhen_df)
  dim(G_full_subset)
  
  #check everyone is lined up
  if(!identical(as.character(subset_indPhen_df$indID), colnames(G_full_subset))){print("Error: not lined up"); break}
  if(!identical(as.character(muts_full$mutname), rownames(G_full_subset))){print("Error: not lined up"); break}
  
  for (i in 1:nrow(G_full_subset)){
    lmgwas_temp <- summary(lm(subset_indPhen_df$phen_temp ~ G_full_subset[i,] + subset_indPhen_df$PC1 + subset_indPhen_df$PC2))
    muts_full$gwas_temp_est[i] <- lmgwas_temp$coeff[2,1]
    muts_full$gwas_temp_P[i] <- lmgwas_temp$coeff[2,4]
    
    lmgwas_sal <- summary(lm(subset_indPhen_df$phen_sal ~ G_full_subset[i,] + subset_indPhen_df$PC1 + subset_indPhen_df$PC2))
    muts_full$gwas_sal_est[i] <- lmgwas_sal$coeff[2,1]
    muts_full$gwas_sal_P[i] <- lmgwas_sal$coeff[2,4]
  }

  #assess TPR of gwas
  muts_full$gwas_sal_sig <- p.adjust(muts_full$gwas_sal_P, method="BH") < 0.05
  (gwas_TPR_sal <- sum(muts_full$gwas_sal_sig & (muts_full$mutSalEffect!=0), na.rm=TRUE)/ # gwas significant and non-zero effect
      sum((muts_full$mutSalEffect!=0), na.rm=TRUE)) #divide by total number of causal mutations
  
  muts_full$gwas_temp_sig <- p.adjust(muts_full$gwas_temp_P, method="BH") < 0.05
  (gwas_TPR_temp <-  sum(muts_full$gwas_temp_sig & (muts_full$mutTempEffect!=0), na.rm=TRUE)/ # gwas significant and non-zero effect
    sum((muts_full$mutTempEffect!=0), na.rm=TRUE)) #divide by total number of causal mutations
  
  (gwas_FDR_sal_neutbase <- sum(muts_full$gwas_sal_sig & muts_full$causal_sal=="neutral")/ #purely neutral loci unaffected by selection
      sum(muts_full$gwas_sal_sig & muts_full$causal_sal!="neutral-linked")) #total number of causal loci and neutral loci outliers
  
  (gwas_FDR_temp_neutbase <- sum(muts_full$gwas_temp_sig & muts_full$causal_temp=="neutral")/ #purely neutral loci unaffected by selection
      sum(muts_full$gwas_temp_sig & muts_full$causal_temp!="neutral-linked")) #total number of causal loci and neutral loci outliers
  
  # GWAS plots
  pdf(paste0(path,seed,"_pdf_7b_GWAS.pdf"), width=7, height=7)
    plot(muts_full$mutTempEffect, muts_full$gwas_temp_est); abline(0,1)
    plot(muts_full$mutSalEffect, muts_full$gwas_sal_est); abline(0,1)
  
    ggplot(muts_full) + geom_density(aes(x=-log10(gwas_temp_P), color=causal_temp)) + ggtheme
    ggplot(muts_full) + geom_density(aes(x=-log10(gwas_sal_P), color=causal_sal)) + ggtheme
    
    ggplot(muts_full) + geom_density(aes(x=abs(gwas_temp_est), color=causal_temp)) + ggtheme
    ggplot(muts_full) + geom_density(aes(x=abs(gwas_sal_est), color=causal_sal)) + ggtheme
    
    ggplot(muts_full[order(muts_full$causal_sal, decreasing=TRUE),]) + geom_point(aes(x=-log10(gwas_sal_P), y=abs(gwas_sal_est), color=causal_sal), alpha=0.5) + ggtheme
    ggplot(muts_full[order(muts_full$causal_temp, decreasing=TRUE),]) + geom_point(aes(x=-log10(gwas_temp_P), y=abs(gwas_temp_est), color=causal_temp), alpha=0.5) + ggtheme
  dev.off()
    
    
  #framework for assessing clinal paradigm ####

  
  criteria_sal <- muts_full$gwas_sal_P < quantile(muts_full$gwas_sal_P,0.05) # top 5% GWAS hits
  sal_gwastop5per_num <- sum(criteria_sal)
  # proportionf of top 5% of GWAS (with structure correction) hits that show clines (without structure correction)
  clinalparadigm_sal_proptop5GWASclines <- sum(criteria_sal & # meets outlier GWAS criteria
        muts_full$cor_sal_sig)/ # AND CLINAL #muts_full$LEA3.2_lfmm2_mlog10P_salenv_sig )/ # This didn't work as well for temp b/c overcorrection for structure
    sal_gwastop5per_num 
  #compare to cor_TPR_sal
  
  # proportion of GWAS outliers (with structure correction) that show clines (without structure correction)
  (clinalparadigm_sal_propsigGWASclines <- sum(muts_full$gwas_sal_sig & # meets outlier GWAS criteria
                                                 muts_full$cor_sal_sig)/ # AND CLINAL #muts_full$LEA3.2_lfmm2_mlog10P_salenv_sig )/ # This didn't work as well for temp b/c overcorrection for structure
    sum(muts_full$gwas_sal_sig))
  #compare to cor_TPR_sal
  
  criteria_temp <- muts_full$gwas_temp_P < quantile(muts_full$gwas_temp_P,0.05)
  temp_gwastop5per_num <- sum(criteria_temp)
  clinalparadigm_temp_proptop5GWASclines <- sum(criteria_temp  & # top 5% of GWAS hits
        muts_full$cor_temp_sig)/ #AND CLINAL #muts_full$LEA3.2_lfmm2_mlog10P_tempenv_sig )/ # This didn't work as well for temp b/c overcorrection for structure
    temp_gwastop5per_num
  #compare to cor_TPR_temp
  
  # proportion of GWAS outliers (with structure correction) that show clines (without structure correction)
  (clinalparadigm_temp_propsigGWASclines <- sum(muts_full$gwas_temp_sig & # meets outlier GWAS criteria
                                                 muts_full$cor_temp_sig)/ # AND CLINAL #muts_full$LEA3.2_lfmm2_mlog10P_salenv_sig )/ # This didn't work as well for temp b/c overcorrection for structure
      sum(muts_full$gwas_temp_sig))
  #compare to cor_TPR_temp
  
  
### AF as a function of environment ####
  
  #str(af_pop)
  #subset_temp_opt
  #subset_sal_opt
  if(!identical(as.character(muts_full$mutname), as.character(rownames(G_full_subset)))){print("Error 1j: mutations not lined up");break()}
  
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
  colnames(af_temp) <- muts_full$mutname
  v <- tapply(af_pop[,1], subset_temp_opt, mean)
  rownames(af_temp) <- as.numeric(dimnames(v)[[1]])
  
  colnames(af_sal) <- muts_full$mutname
  v <- tapply(af_pop[,1], subset_sal_opt, mean)
  rownames(af_sal) <- as.numeric(dimnames(v)[[1]])

  muts_full$af_cor_temp_pooled <- as.numeric(cor(af_temp, temp_levels, method = "kendall"))
  muts_full$af_cor_sal_pooled <- as.numeric(cor(af_sal, sal_levels, method = "kendall"))

  pdf(paste0(path,seed,"_pdf_7_afcors.pdf"), width=9, height=8)
  
  ggplot() + 
    geom_point(data=muts_full[muts_full$causal_temp=="neutral-linked",], aes(x=af_cor_temp, y = af_cor_temp_pooled), color=adjustcolor("goldenrod", 0.1)) +
    geom_point(data=muts_full[muts_full$causal_temp=="neutral",], aes(x=af_cor_temp, y = af_cor_temp_pooled), color=adjustcolor("grey", 0.1)) +
    geom_point(data = muts_full[muts_full$causal_temp=="causal",], 
               aes(x=af_cor_temp, y = af_cor_temp_pooled,  size=Va_temp_prop, color=Va_temp_prop), alpha=0.8, shape=shape_causal) +
    scale_colour_viridis(option=cs, begin = begin_cs, end=end_cs, limits=c(0,1)) +  ggtheme +  xlab("Cor(af, temp)") + 
    ylab("Cor(af, temp) pooled by temp.") + ggtitle(paste0(plotmain, " mutations: temp")) + ylim(-1,1) + xlim(-1, 1) + geom_abline(data=NULL, intercept=0, slope=1) + labs(color="temp VA prop.", size="temp VA prop.")
    
  ggplot() + 
    geom_point(data=muts_full[muts_full$causal_sal=="neutral-linked",], aes(x=af_cor_sal, y = af_cor_sal_pooled), color=adjustcolor("goldenrod", 0.1)) +
    geom_point(data=muts_full[muts_full$causal_sal=="neutral",], aes(x=af_cor_sal, y = af_cor_sal_pooled), color=adjustcolor("grey", 0.1)) +
    geom_point(data = muts_full[muts_full$causal_sal=="causal",], 
               aes(x=af_cor_sal, y = af_cor_sal_pooled,  size=Va_sal_prop, color=Va_sal_prop), alpha=0.8, shape=shape_causal) +
    scale_colour_viridis(option=cs, begin = begin_cs, end=end_cs, limits=c(0,1)) +  ggtheme +  xlab("Cor(af, Env2)") + 
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

  #col <- two.colors(num_causal_postfilter, start="grey", middle="goldenrod", end="orange", alpha=0.7)
  col <- turbo(num_causal_postfilter, begin=0.8, end=0)
  #col <- viridis(num_causal_postfilter, begin=1, end=0)
  #col <- mako(num_causal_postfilter, begin=1, end=0.2)
  
  muts_full$color_af.temp.cline <- muts_full$color_af.sal.cline <- NA
  muts_full1 <- muts_full[order(abs(muts_full$cor_temp), decreasing=TRUE),]
  muts_full1$color_af.temp.cline[1:length(col)] <- as.character(col)
  muts_full <- muts_full1[order(muts_full1$VCFrow),]
  
  if(!identical(as.character(muts_full$mutname), as.character(rownames(G_full_subset)))){print("Error 1k: mutations not lined up");break()}
  
  
  muts_full1 <- muts_full[order(abs(muts_full$cor_sal), decreasing=TRUE),]
  muts_full1$color_af.sal.cline[1:length(col)] <- as.character(col)
  muts_full <- muts_full1[order(muts_full1$VCFrow),]
  
  if(!identical(as.character(muts_full$mutname), as.character(rownames(G_full_subset)))){print("Error 1k: mutations not lined up");break()}
  
  
  #Temp plot
  plot(pop_df$opt1, rep(0, length(pop_df$opt1)), ylim=c(0,1), xlab="Deme Temperature", 
       ylab="Allele frequency", col=rgb(0,0,0,0), bty="n", xlim=c(-1, 1), las=1)
  
  #str(af_temp)
  for (i in which(muts_full$causal_temp=="causal")){
    cor_mut <- muts_full$cor_temp[i]
    if(abs(cor_mut)<0.2){
      lwd=5
      lty=1
      alpha1 = 0.5
    }
    if(abs(cor_mut)>0.2 & abs(cor_mut)<0.5){
      lwd=3
      lty=1
      alpha1 = 0.5
    }
    if(abs(cor_mut)>0.5){
      lwd=2
      lty=2
      alpha1 = 0.5
    }
    lines(temp_levels, af_temp[,i], 
          col=adjustcolor(muts_full$color_af.temp.cline[i], alpha1), lwd=lwd, lty=lty)
  }  

  plot(pop_df$opt0, rep(0, length(pop_df$opt0)), ylim=c(0,1), xlab="Deme Env2", 
       ylab="Allele frequency", col=rgb(0,0,0,0), bty="n", xlim=c(-1, 1), las=1)  
  for (i in which(muts_full$causal_sal=="causal")){
    cor_mut <- muts_full$cor_sal[i]
    if(abs(cor_mut)<0.2){
      lwd=5
      lty=1
      alpha1 = 0.5
    }
    if(abs(cor_mut)>0.2 & abs(cor_mut)<0.5){
      lwd=3
      lty=1
      alpha1 = 0.5
    }
    if(abs(cor_mut)>0.5){
      lwd=2
      lty=2
      alpha1 = 0.5
    }
    lines(sal_levels, af_sal[,i], col=adjustcolor(muts_full$color_af.sal.cline[i], 
                                                  alpha1), lwd=lwd, lty=lty)
  }  

  ## Make the legends ##
  # Create two panels side by side
  layout(t(1:2), widths=c(5,1))
  # Set margins and turn all axis labels horizontally (with `las=1`)
  par(mar=c(0.5, 0.5,4,0.5), oma=c(0, 0, 0,4), las=1)
  
  plot(pop_df$opt1, rep(0, length(pop_df$opt1)), ylim=c(0,1), xlab="Temperature", 
       ylab="Allele frequency", 
       col=rgb(0,0,0,0), bty="n", xlim=c(-1, 1), main=plotmain, cex.main=0.5)
  
  legend(-1,0.5, c("linear model", "1:1 line") , lty=c(1,2), bty="n")
  legend(0, 0.5, c("abs(cor) < 0.2", "0.2 < abs(cor) > 0.5", 
                     "abs(cor) > 0.5"), lwd=c(5, 3, 2), lty=c(1,1,2), bty="n")

  # Draw the color legend
  image(1, seq(0.0,1, length.out=100), t(seq_along(seq(0,1, length.out=100))), 
        col=rev(col), axes=FALSE, las=1,
        main="abs[Cor(p,env)]", cex.main=0.7)
  axis(4)
  
  dev.off()
  #legend(0,1, seq(0.05,0.95, length.out=10), fill=turbo(10, begin = 0, end=0.8), 
  #       cex=0.7, bty="n", title="abs(Cor (p, env))") 
  
  
### genotype heatmaps ####
  
pdf(paste0(path,seed,"_pdf_heatmaps.pdf"), width=8, height=8)

par(cex.main=0.5)

# Make sure in correct order
if(!identical(subset_indPhen_df$temp_opt , subset_indPhen_df$temp_opt[order(subset_indPhen_df$temp_opt)])){print("Error optimums not identical");break}

G_heatmap <- G_full_subset[which(muts_full$causal),]

a<-heatmap(t(G_heatmap), Rowv = NA,  
           main=paste0("\n",plotmain,"All Causal SNPs"),cexCol = 0.3,
           #Colv = NA,
           labRow = round(subset_indPhen_df$temp_opt,2),
           ylab="Temperature (South <---------> North)", cex.main=0.5, scale="none")
a

# low salinity sites
heatmap(t(G_heatmap[a$colInd,subset_indPhen_df$sal_opt==-1]), Rowv = NA,  
        main=paste0("\n",plotmain, "All causal SNPs, Low Env2 = -1 Sites"),cexCol = 0.3,
        Colv = NA,
        labRow = round(subset_indPhen_df$temp_opt[subset_indPhen_df$sal_opt==-1],2),
        ylab="Temperature (South <---------> North)", cex.main=0.5, scale="none")

heatmap(t(G_heatmap[a$colInd,subset_indPhen_df$sal_opt==1]), Rowv = NA,  
        main=paste0("\n",plotmain, "All causal SNPs, High Env2 = 1 Sites"),cexCol = 0.3,
        Colv = NA,
        labRow = round(subset_indPhen_df$temp_opt[subset_indPhen_df$sal_opt==1],2),
        ylab="Temperature (South <---------> North)", cex.main=0.5, scale="none")

heatmap(t(G_heatmap[a$colInd,abs(subset_indPhen_df$sal_opt)<0.15]), Rowv = NA,  
        main=paste0("\n",plotmain, "All causal SNPs, Intermediate Env2 ~ 0 Sites"),cexCol = 0.3,
        Colv = NA,
        labRow = round(subset_indPhen_df$temp_opt[abs(subset_indPhen_df$sal_opt)<0.15],2),
        ylab="Temperature (South <---------> North)", cex.main=0.5, scale="none")

heatmap(t(G_heatmap[a$colInd,subset_indPhen_df$x==1]), Rowv = NA,  
        main=paste0("\n",plotmain, "All causal SNPs, x = 1 Sites"),cexCol = 0.3,
        Colv = NA,
        labRow = round(subset_indPhen_df$temp_opt[subset_indPhen_df$x==1],2),
        ylab="Temperature (South <---------> North)", cex.main=0.5, scale="none")

heatmap(t(G_heatmap[a$colInd,subset_indPhen_df$x==10]), Rowv = NA,  
        main=paste0("\n",plotmain, "All causal SNPs, x = 10 Sites"),cexCol = 0.3,
        Colv = NA,
        labRow = round(subset_indPhen_df$temp_opt[subset_indPhen_df$x==10],2),
        ylab="Temperature (South <---------> North)", cex.main=0.5, scale="none")


heatmap(t(G_heatmap[a$colInd,subset_indPhen_df$x==5]), Rowv = NA,  
        main=paste0("\n",plotmain, "All causal SNPs, x = 5 Sites"),cexCol = 0.3,
        Colv = NA,
        labRow = round(subset_indPhen_df$temp_opt[subset_indPhen_df$x==5],2),
        ylab="Temperature (South <---------> North)", cex.main=0.5, scale="none")



dev.off()

### Write outputs ####

write.table(muts_full, paste0(path,seed,"_Rout_muts_full.txt"))
write.table(subset_indPhen_df, paste0(path,seed,"_Rout_ind_subset.txt"))
write.table(af_pop, paste0(path,seed,"_Rout_af_pop.txt"))
write.table(af_sal, paste0(path,seed,"_Rout_af_sal.txt"))
write.table(af_temp, paste0(path,seed,"_Rout_af_temp.txt"))
            
# cor_VA_temp
# cor_VA_sal
# cor_TPR_temp
# cor_FDR_temp
# cor_TPR_sal
# cor_FDR_sal

out_full <- data.frame(seed=as.character(seed),
                       n_samp_tot=nrow(subset_indPhen_df), 
                       n_samp_per_pop = n_per_pop,
                       sd_fitness_among_inds,
                       sd_fitness_among_pops,
                       final_LA = LA_df$local_adapt[nrow(LA_df)],
                       K, Bonf_alpha,
                       numCausalLowMAFsample,
                       all_corr_phen_temp,
                       subsamp_corr_phen_temp,
                       all_corr_phen_sal,
                       subsamp_corr_phen_sal,
                       num_causal_prefilter,
                       num_causal_postfilter,
                       num_non_causal = sum(!muts_full$causal),
                       num_neut_prefilter = num_neut_prefilter, #total number non-causal loci
                       num_neut_postfilter = num_neut_prefilter, #total number non-causal loci after filtering
                       num_neut_neutralgenome = sum(muts_full$causal_temp=="neutral"), #loci in 2nd half of genome
                       num_causal_temp = sum(muts_full$causal_temp=="causal"),
                       num_causal_sal = sum(muts_full$causal_sal=="causal"),
                       num_multiallelic,
                       
                       meanFst,
                       va_temp_total,
                       va_sal_total,
                       
                       Va_temp_sample,
                       Va_sal_sample,
                       nSNPs, #` total number of SNPs in analysis
                      
                      median_causal_temp_cor, # median correlation between AF and Temp for causal loci
                       median_causal_sal_cor, 
                       median_neut_temp_cor, 
                       median_neut_sal_cor,           
                       
                # CORRELATION OUTPUTS       
                       cor_VA_temp_prop,
                       cor_VA_sal_prop,
                       cor_TPR_temp,
                       cor_TPR_sal,
                       
                       cor_FDR_allSNPs_temp,
                       cor_FDR_neutSNPs_temp,
                       cor_FDR_allSNPs_sal,
                       cor_FDR_neutSNPs_sal,
                      num_causal_sig_temp_corr, ## number of causal loci that have significant Kendall's correlations with temperature after Bonferroni correction
                      num_causal_sig_sal_corr, ## number of causal loci that have significant Kendall's correlations with salinity after Bonferroni correction
                      num_notCausal_sig_temp_corr, # number of non-causal (neutral and neutral-linked) loci that have significant Kendall's correlations with temperature after Bonferroni correction
                      num_notCausal_sig_sal_corr, # number of non-causal (neutral and neutral-linked) loci that have significant Kendall's correlations with salinity after Bonferroni correction
                      num_neut_sig_temp_corr, # number of neutral (unlinked to causal) loci that have significant Kendall's correlations with temperature after Bonferroni correction
                      num_neut_sig_sal_corr,# number of neutral (unlinked to causal) loci that have significant Kendall's correlations with salinity after Bonferroni correction
                             
                       cor_AUCPR_temp_allSNPs,
                       cor_AUCPR_temp_neutSNPs,
                       cor_AUCPR_sal_allSNPs,
                       cor_AUCPR_sal_neutSNPs,
                       cor_af_temp_noutliers, # number of outliers for cor(af,temp) after Bonferroni correction
                       cor_af_sal_noutliers, # number of outliers for cor(af,salinity) after Bonferroni correction
                       cor_FPR_temp_neutSNPs, #` false positive rate in cor(af,temp) after Bonferroni correction, based on loci unaffected by selection
                       cor_FPR_sal_neutSNPs, #` false positive rate in cor(af,sal) after Bonferroni correction, based on loci unaffected by selection
                      
                # LFMM OUTPUTS 
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
                      LEA3.2_lfmm2_mlog10P_tempenv_noutliers, #` number of outliers for the lfmm temp model (qvalue <0.05)
                      LEA3.2_lfmm2_mlog10P_salenv_noutliers, #` number of outliers for the lfmm salinity model (qvalue <0.05)
                      LEA3.2_lfmm2_num_causal_sig_temp, #` number of causal loci on the temp trait, significant in the lfmm temp model (qvalue <0.05)
                      LEA3.2_lfmm2_num_neut_sig_temp, #` number of neutral loci false positives (only neutral loci not affected by selection), significant in the lfmm temp model (qvalue <0.05)
                      LEA3.2_lfmm2_num_causal_sig_sal, #` number of causal loci on the salinity trait, significant in the lfmm salinity model (qvalue <0.05)
                      LEA3.2_lfmm2_num_neut_sig_sal,
                      LEA3.2_lfmm2_FPR_neutSNPs_temp, #` false positive rate of lfmm temperature model
                      LEA3.2_lfmm2_FPR_neutSNPs_sal, #` false positive rate of lfmm salinity model
                
                #RDA outputs
                      RDA1_propvar, #proportion of variance explained by 1st RDA axis
                      RDA2_propvar, #proportion of variance explained by 2nd RDA axis
                      RDA1_propvar_corr, #proportion of variance explained by 1st RDA axis with structure correction
                      RDA2_propvar_corr, #proportion of variance explained by 2nd RDA axis with structure correction
                      
                      RDA1_temp_cor, #` output of `summary(rdaout)$biplot[2,1]`, which is the correlation between RDA1 and the temperature environmental variable
                      RDA1_sal_cor, #`  output of `summary(rdaout)$biplot[1,1]`, which is the correlation between RDA1 and the salinity environmental variable
                      RDA2_temp_cor, #` output of `summary(rdaout)$biplot[2,2]`, which is the correlation between RDA2 and the temperature environmental variable
                      RDA2_sal_cor, #` output of `summary(rdaout)$biplot[1,2]`, which is the correlation between RDA2 and the salinity environmental variable
                      
                # RDA OUTLIER OUTPUTS
                      RDA_Va_temp_prop,
                      RDA_Va_temp_prop_corr, #prop of VA explained by outliers with structure correction
                      RDA_Va_sal_prop,
                      RDA_Va_sal_prop_corr, #prop of VA explained by outliers with structure correction
                      RDA_TPR, 
                      RDA_TPR_corr, #with structure correction
                      RDA_FDR_allSNPs,
                      RDA_FDR_allSNPs_corr,#with structure correction
                      num_RDA_sig_causal,
                      num_RDA_sig_neutral, #only for neutral half of genome
                      num_RDA_sig_causal_corr, #with structure correction
                      num_RDA_sig_neutral_corr, #with structure correction
                      RDA_FDR_neutSNPs,
                      RDA_FDR_neutSNPs_corr,  #with structure correction
                      RDA_AUCPR_allSNPs,
                      RDA_AUCPR_neutSNPs,
                      RDA_AUCPR_neutSNPs_corr, #with structure correction
                      RDA_FPR_neutSNPs, 
                      RDA_FPR_neutSNPs_corr, #with structure correction
          
                # RDA PREDICTION OUTPUTS     
                      RDA_RDAmutpred_cor_tempEffect, #kendall's correlation between the predicted temperature effect from RDA and the true mutation effect on temperature
                      RDA_RDAmutpred_cor_salEffect, # kendall's correlation between the predicted salinity effect from RDA and the true mutation effect on salinity
                      RDA_absRDAmutpred_cor_tempVa, # kendall's correlation between the abs(predicted temperature effect from RDA) and the true mutation Va on temperature
                      RDA_absRDAmutpred_cor_salVa, #kendall's correlation between the abs(predicted salinity effect from RDA) and the true mutation Va on salinity
                    
                      RDA_RDAmutpred_cor_tempEffect_structcorr,
                      RDA_RDAmutpred_cor_salEffect_structcorr,
                      RDA_absRDAmutpred_cor_tempVa_structcorr,
                      RDA_absRDAmutpred_cor_salVa_structcorr,
                      
                      RDA_cor_RDA20000temppredict_tempPhen, 
                      RDA_cor_RDA20000salpredict_salPhen,
                      RDA_cor_RDA20000temppredict_tempPhen_structcorr,
                      RDA_cor_RDA20000salpredict_salPhen_structcorr,
                

            # STRUCUTRE CORRELATION OUTPUTS
                      cor_PC1_temp, 
                      cor_PC1_sal, 
                      cor_PC2_temp, 
                      cor_PC2_sal, 
                      cor_LFMMU1_temp, 
                      cor_LFMMU1_sal, 
                      cor_LFMMU2_temp, 
                      cor_LFMMU2_sal,
                      cor_PC1_LFMMU1_temp, 
                      cor_PC1_LFMMU1_sal, 
                      cor_PC2_LFMMU1_temp, 
                      cor_PC2_LFMMU1_sal,
            
            #GWAS outputs
            gwas_TPR_sal,
            gwas_TPR_temp,
            gwas_FDR_sal_neutbase,
            gwas_FDR_temp_neutbase,
            clinalparadigm_sal_proptop5GWASclines,
            clinalparadigm_temp_proptop5GWASclines,
            clinalparadigm_sal_propsigGWASclines,
            clinalparadigm_temp_propsigGWASclines
)

write.table(out_full,paste0(path,seed,"_Rout_simSummary.txt"), row.names = FALSE)


  
