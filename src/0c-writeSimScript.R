setwd("/Users/lotterhos/Documents/GitHub/MVP-NonClinalAF")

final <- read.table("src/0b-final_params.txt", header=TRUE)
head(final)

#final[1:15,]
#levels(factor(final$level))
levels(factor(final$demog_name))
levels(factor(final$arch))

head(final[grep("oliogenic_2-trait-pleiotropy-equal-S", final$arch),], 15)

head(final[grep("oliogenic_2-trait-pleiotropy-equal-S__Est-Clines_N-equal_m-constant", final$level),],1)
which(final$seed==1231142) # done check

# To check
head(final[grep("polygenic_2-trait-pleiotropy-equal-S__Est-Clines_N-equal_m-constant", final$level),],1)
which(final$seed==1231217)

head(final[grep("oliogenic_2-trait-pleiotropy-equal-S__Est-Clines_N-equal_m_breaks", final$level),],1)
which(final$seed==1231141)

head(final[grep("oliogenic_2-trait-pleiotropy-equal-S__SS-Clines_N-equal_m_breaks", final$level),],1)
which(final$seed==1231146)

head(final[grep("oliogenic_2-trait-pleiotropy-unequal-S__Est-Clines_N-variable_m-variable", final$level),],1)
which(final$seed==1231158)

head(final[grep("oliogenic_2-trait-pleiotropy-unequal-S__Est-Clines_N-equal_m-constant", final$level),],1)
which(final$seed==1231157)


i=53
final[i,]
final$seed[i]
final$level[i]
outfile=paste0(final$seed[i],"_outfile.txt")
outpath="sim_outputs/"
outpathfile=paste0(outpath, outfile)

final$level[i]

paste0("cat > sim_outputs/", final$seed[i],"__", final$level[i], ".txt")

# assumes path is set to MVP-NonClinalAF folder
cat("slim -d MY_SEED1=",final$seed[i],
    " -d \"MY_RESULTS_PATH1='sim_outputs/'\"",
    " -d MIG_x1=", final$MIG_x[i],
    " -d MIG_y1=",final$MIG_y[i],
    " -d \"demog1='",final$demog[i],
    "'\" -d \"xcline1='", final$xcline[i],
    "'\" -d \"ycline1='" , final$ycline[i],
    "'\" -d METAPOP_SIDE_x1=",final$METAPOP_SIDE_x[i] ,
    " -d METAPOP_SIDE_y1=", final$METAPOP_SIDE_y[i],
    " -d Nequal1=",final$Nequal[i],
    " -d isVariableM1=", final$isVariableM[i],
    " -d MIG_breaks1=" ,final$MIG_breaks[i],
    " -d MU_base1=", final$MU_base[i],
    " -d MU_QTL_proportion1=", final$MU_QTL_proportion[i],
    " -d SIGMA_QTN_1a=", final$SIGMA_QTN_1[i],
    " -d SIGMA_QTN_2a=",final$SIGMA_QTN_2[i],
    " -d SIGMA_K_1a=", final$SIGMA_K_1[i],
    " -d SIGMA_K_2a=", final$SIGMA_K_2[i],
    " -d N_traits1=", final$N_traits[i],
    " -d ispleiotropy1=", final$ispleiotropy[i],
    " src/a-PleiotropyDemog.slim",
    " > ", outpathfile," 2> ", outpathfile, ".error.txt", sep="")

cat("python3 src/b-process_trees.py -s ",
    final$seed[i],
    " -r 1e-06",
    " -mu ", final$MU_base[i]*10,
    " -N 1000 > ", outpath, final$seed[i],
    "_pytree.out.txt 2> ",
    outpath, final$seed[i], "_pytree.error.txt", sep="")
