
setwd("/work/lotterhos/MVP-NonClinalAF")

# The loop takes about 10 hours to run

### load params ####
  
  #args = commandArgs(trailingOnly=TRUE)
  
  #seed = 1231094
  #seed = 1231098
  # seed=1231139
  #seed = 1231144
  #seed=1232947
  #seed = 1231102
  # seed = 1231109
  #path = "sim_output_20220428/"
  #runID=20220428
  #seed = args[1]
  #path = args[2]
  #runID = args[3]

### Read in Sims
allsims <- load(paste0("src/0b-final_params-",runID,".RData"))
allsims<- final

which(allsims$seed==1231109)

for (seed in allsims$seed[16:nrow(allsims)]){
  #seed = 1231109
  thissim <- allsims[grep(seed, allsims$seed),]
  (plotmain <- paste(thissim$level, seed, sep="\n"))
  print(plotmain)
  
  ### Read in DAta ####
  info <- read.table(paste0(path,seed,"_info.txt"), header=TRUE, 
                     colClasses = c("character", 
                                    rep("numeric",15),
                                    "character",
                                    rep("numeric",6)))
  info
  
  out <- read.table(paste0(path,seed,"_Rout_simSummary.txt"), header=TRUE)
  out$final_LA
  
  G <- read.table(paste0(path,seed,"_Rout_Gmat_sample.txt"), header=TRUE)
  muts <- read.table(paste0(path,seed,"_Rout_muts_full.txt"), header=TRUE)
  ind <- read.table(paste0(path,seed,"_Rout_ind_subset.txt"), header=TRUE)
  pop_df <- read.table(paste0(path,seed,"_popInfo.txt"), header=TRUE, 
                       colClasses = c("character", rep("numeric",6))) 
  
  head(G[,1:100])  
  
  if (nrow(G)!=nrow(muts)){
    print(paste("Error in seed", seed, "number of mutations in data frame does not match G matrix"))
  }
  #dim(G)
  #dim(muts)
  #dim(ind)
  
  # just a check that multiplication works how I think
  # a <- matrix(sample(c(0,1), 100, replace=TRUE), nrow=50)    
  # a
  # x <- seq(1:50)  
  # a*x
  
  
  
  ### Function to calculate the fitness of an individual given the stabilizing selection function ####
  fitness_var_temp <- info$SIGMA_K_2
  fitness_var_sal <- info$SIGMA_K_1
  fitness_cov <- info$SIGMA_K_Cov
  N_traits <- info$Ntraits
  
  get_fitness <- function( temp_phen, temp_opt, sal_phen=NA, sal_opt=NA){
    # For one trait, temp_phen can be a vector and temp_opt can be a corresponding vector
    # for two traits, the phenotypes can be vectors of multiple individuals, but the optimum is a single value
    
    if (N_traits==2){   # 2 traits
      phens <- cbind(temp_phen, sal_phen)
      opts <- c(temp_opt, sal_opt)
      fitness_varcov <- matrix(c(fitness_var_temp, fitness_cov, fitness_cov, fitness_var_sal), nrow=2);
      fitness_norm <- dmvnorm(c(0.0, 0.0), c(0.0, 0.0), fitness_varcov);
      fits <- dmvnorm(phens, opts, fitness_varcov) / fitness_norm;
    }
    
    if (N_traits==1){ # 1 trait
      fitness_norm = dnorm(0.0, 0.0, fitness_var_temp);
      fits = dnorm(temp_phen, temp_opt, fitness_var_temp)/fitness_norm;
      
    }
    return(fits)
    
  }
  
  
  ### Function to calculate the amount of local adaptation given individuals, phentypes, and locations ####
  amount_LA <- function(temp_phen, sal_phen, popIDs){
    # temp_phen and sal_phen are vectors of phenotypes
    
    fitness_matrix_pop = matrix(rep(0.0, 100*100), nrow=100, ncol=100);
    # the row in this matrix is the source genotype
    # the column in this matrix is the transplant population
    
    # for the first transplant population, get the mean genotype from all the other source populations
    # calculate their fitness at the transplant location
    # fill in the column of the matrix with the genotype fitnesses
    
    for (i in 1:100){ # hard coded for 100 populations
      # i is the transplat site
      temp_phen_mean <- tapply(temp_phen, popIDs, mean)
      sal_phen_mean <- tapply(sal_phen, popIDs, mean)
      #pop_df$opt0[i] # optimum salinity phenotype at transplant site
      #pop_df$opt1[i] # optimum temp phenotype at transplant site
      
      if (N_traits ==1){
        fitness_matrix_pop[,i] <- get_fitness(temp_phen_mean, pop_df$opt1[i])
      }
      
      if (N_traits ==2){
        fitness_matrix_pop[,i] <- get_fitness(temp_phen_mean, pop_df$opt1[i], sal_phen_mean , pop_df$opt0[i])
      }
    } # end loop through transplant sites
    
    diagonals <- rep(NA, 100)
    for (i in 1:100){ # calculate LA
      diagonals[i] = fitness_matrix_pop[i,i];
    } # end for loop
    
    sympatric = mean(diagonals);
    allopatric = (sum(fitness_matrix_pop) - sum(diagonals)) / (length(fitness_matrix_pop) - length(diagonals));
    local_adapt = sympatric - allopatric;
    
    return(local_adapt)
  } # end function
  
  ### Total LA for metapopulation ####
  out$final_LA
  
  ### Calculate amount of LA for sampled individuals using their phenotypes from the simulation with all loci ####
  (LAcalc_1000ind_allmuts <- amount_LA(ind$phen_temp, ind$phen_sal, ind$subpopID))
  
  ### Calculate amount of LA for sample (MAF > 0.01) and compare to amount calculated with all mutations in simulation ####
  
  effectG_temp <- G*muts$mutTempEffect
  effectG_sal <- G*muts$mutSalEffect
  
  par(mfrow=c(2,1), mar=c(4,4,1,1))
  ind$phen_temp_mutsSample <- colSums(effectG_temp, na.rm=TRUE)
  #plot(ind$phen_temp, ind$phen_temp_mutsSample)  # sanity check
  #abline(0,1, col="orange")
  
  ind$phen_sal_mutsSample <-colSums(effectG_sal, na.rm=TRUE)
  #plot(ind$phen_sal, ind$phen_sal_mutsSample)  # sanity check
  #abline(0,1, col="orange")
  
  LAcalc_1000ind_NoMAFmuts <- amount_LA(ind$phen_temp_mutsSample, ind$phen_sal_mutsSample, ind$subpopID)
  
  
  ### Calculate amount of LA for clinal alleles in both environments and compare to amount calculated with all mutations in simulation ####
  alpha <- 0.05/nrow(G)
  keep <- (muts$af_cor_temp_P < alpha) | (muts$af_cor_sal_P < alpha)
  
  if (sum(keep) > 0){
    ind$phen_temp_mutsClinal <- colSums(effectG_temp[which(keep),], na.rm=TRUE)
    ind$phen_sal_mutsClinal <- colSums(effectG_sal[which(keep),], na.rm=TRUE)    
    LAcalc_1000ind_ClinalMuts <- amount_LA(ind$phen_temp_mutsClinal , ind$phen_sal_mutsClinal , ind$subpopID)
    
    ### Calculate amount of LA for random set of loci same number as was used in the last one ####
    xs <- sample(1:nrow(G), sum(keep))
    ind$phen_temp_mutsRandom <- colSums(effectG_temp[xs,], na.rm=TRUE)
    ind$phen_sal_mutsRandom <- colSums(effectG_sal[xs,], na.rm=TRUE)
    LAcalc_1000ind_nClinal <- sum(keep)
    LAcalc_1000ind_RandomMuts_nClinal <- amount_LA(ind$phen_temp_mutsRandom , ind$phen_sal_mutsRandom , ind$subpopID)
  }else{
    LAcalc_1000ind_nClinal <-0
    LAcalc_1000ind_ClinalMuts <- 0
    LAcalc_1000ind_RandomMuts_nClinal <- 0
  }
  
  
  
  ### Calculate amount of LA for non-clinal alleles in both environments ####
  ind$phen_temp_mutsNonClinal <- colSums(effectG_temp[which(!keep),], na.rm=TRUE)
  ind$phen_sal_mutsNonClinal <- colSums(effectG_sal[which(!keep),], na.rm=TRUE)
  LAcalc_1000ind_nNonClinal <- sum(!keep)
  LAcalc_1000ind_NonClinalMuts <- amount_LA(ind$phen_temp_mutsNonClinal , ind$phen_sal_mutsNonClinal , ind$subpopID)
  
  ### Calculate amount of LA for random set of loci same number as was used in the last one ####
  xs <- sample(1:nrow(G), sum(!keep))
  ind$phen_temp_mutsRandom <- colSums(effectG_temp[xs,], na.rm=TRUE)
  ind$phen_sal_mutsRandom <- colSums(effectG_sal[xs,], na.rm=TRUE)
  LAcalc_1000ind_RandomMuts_nNonClinal <- amount_LA(ind$phen_temp_mutsRandom , ind$phen_sal_mutsRandom , ind$subpopID)
  
  ### Calculate amount of LA for DETECTABLE Clinal alleles with LFMM    ####
  keep_LFMM <- (muts$LEA3.2_lfmm2_mlog10P_tempenv_sig) | (muts$LEA3.2_lfmm2_mlog10P_salenv_sig)
  
  if (sum(keep_LFMM)>0){
    ind$phen_temp_mutsClinal_LFMM <- colSums(effectG_temp[which(keep_LFMM),], na.rm=TRUE)
    ind$phen_sal_mutsClinal_LFMM <- colSums(effectG_sal[which(keep_LFMM),], na.rm=TRUE)   
    LAcalc_1000ind_ClinalMutsLFMM <- amount_LA(ind$phen_temp_mutsClinal_LFMM , ind$phen_sal_mutsClinal_LFMM, ind$subpopID)    
  }else{
    LAcalc_1000ind_ClinalMutsLFMM <- 0
  }
  
  
  ### Outputs ####
  LAcalc_out <- data.frame(LAcalc_SLiM = out$final_LA,
                           LAcalc_1000ind_allmuts,
                           LAcalc_1000ind_NoMAFmuts,
                           LAcalc_1000ind_ClinalMuts,
                           LAcalc_1000ind_ClinalMutsLFMM,
                           LAcalc_1000ind_NonClinalMuts ,
                           LAcalc_1000ind_RandomMuts_nClinal,
                           LAcalc_1000ind_RandomMuts_nNonClinal
  )        
 # par(mar=c(12,4,2,1), mfrow=c(1,1))
 #  boxplot(LAcalc_out, las=2, ylab="LA estimate", bty="n", border=c(rep("blue",4),"red","blue","grey", "grey"),
 #        ylim=c(0,0.6), main=plotmain, cex.main=0.5)    
  
  causal <- muts$mutSalEffect != 0 | muts$mutTempEffect != 0
  LAcalc_out <- data.frame( LAcalc_out,
                            LAcalc_1000ind_nClinal = sum(keep & causal, na.rm=TRUE),
                            LAcalc_1000ind_nNonClinal = sum(!keep & causal, na.rm=TRUE) , seed)
  LAcalc_out <- round(LAcalc_out,3)
  
  if (seed == allsims$seed[1]){
    write.table(LAcalc_out,"LAcalc.txt", row.names = FALSE, col.names=TRUE)
  }else{
    write.table(LAcalc_out,"LAcalc.txt", row.names = FALSE, col.names=FALSE, append = TRUE)
  }
  
}

