seed = 1231096 # 1231094 to 1231098 for Estuary; 1231099 to 1231103 for SS

path= "sim_output_20220201/"

for (seed in c(1231094, 1231104)){
  pop_df <- read.table(paste0(path,seed,"_popInfo.txt"), header=TRUE, 
                       colClasses = c("character", rep("numeric",6))) 
  
  pdf(paste0("graphs_mig/",seed,"_pdf_1pop_opt.pdf"), width=4, height=4)  
    ggplot(pop_df) + ggtheme + geom_point(aes(x=x, y=y,size=opt0), color="grey20") + geom_point(aes(x=x, y=y, color=opt1), size=2.5) + scale_colour_gradient2(high=rgb(1,0.4,0.2), low="cornflowerblue", mid=rgb(0.8,0.8,0.7), name="Temp") +  labs(size="Env2") + ggtitle(plotmain) + theme(plot.title = element_text(size = 7))
    ggplot(pop_df) + ggtheme + geom_point(aes(x=x, y=y,size=opt0), color="grey20") + geom_point(aes(x=x, y=y, color=opt1), size=2.5) + scale_colour_gradient2(high=rgb(1,0.4,0.2), low="cornflowerblue", mid=rgb(0.8,0.8,0.7), name="Temp")  + 
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + theme(legend.position = "none")
  dev.off()
}

for (seed in 1231094:1231103){

pop_df <- read.table(paste0(path,seed,"_popInfo.txt"), header=TRUE, 
                     colClasses = c("character", rep("numeric",6))) 

allsims <- load(paste0("src/0b-final_params-",runID,".RData"))
allsims<- final
thissim <- allsims[grep(seed, allsims$seed),]
(plotmain <- paste(thissim$level, seed, sep="\n"))
mig <- read.table(paste0(path,seed,"_popInfo_m.txt"), header=TRUE)


mig_thick <- rep(NA, length(mig$m))
mig_thick[mig$m<0.03] <- 0.1
mig_thick[mig$m>=0.03 & mig$m<0.07] <- 1.5
mig_thick[mig$m>=0.07 & mig$m<0.15] <- 3
mig_thick[mig$m>=0.15 & mig$m<0.3] <- 5
mig_thick[mig$m>=0.3] <- 7

par(mfrow=c(1,1), mar=c(4,4,3,1))
pdf(paste0("graphs_mig/",seed,"_pdf_1pop_mig.pdf"), width=8, height=8)  
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
}