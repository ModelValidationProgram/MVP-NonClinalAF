
# Spatial area

# SW 49.923502, -123.191617
# NE 53.407794, -114.808504
library(raster)
library(gplots)
setwd("~/Documents/GitHub/MVP-NonClinalAF/multipheno_multienvi/bioclim")

minlat <- 49.923
maxlat <- 53.408
gridsize <- 360
minlong <- -123.192
maxlong <- -114.809
lats <- round(seq(minlat, maxlat, by= (maxlat-minlat)/(gridsize-1)),3)
lons <- round(seq(minlong, maxlong, by= (maxlong-minlong)/(gridsize-1)),3)

lat_grid <- rep(lats, gridsize)
long_grid <- rep(lons, each=gridsize)

## Altitude
a <- getData("worldclim",var="alt",res=2.5)

# WorldClim
r <- getData("worldclim",var="bio",res=2.5)
#head(r)

bios <- data.frame(BIO=paste0("BIO",1:19), 
                   ABBRV= c("MAT", "MDR", "ISO",
                                                   "TSsd", "MaxTWM", "MTCM",
                                                   "TAR", "MTWetQ", "MTDQ",
                                                   "MTWarmQ", "MTCQ", "MAP",
                                                   "PWM", "PDM", "PSsd",
                                                   "PWwetQ", "PDQ", "PWarmQ", "PCQ"),
                   DESC =c("Annual Mean Temperature",
"Mean Diurnal Range (Mean of monthly (max temp - min temp))",
"Isothermality (BIO2/BIO7) (×100)",
"Temperature Seasonality (standard deviation ×100)",
"Max Temperature of Warmest Month",
"Min Temperature of Coldest Month",
"Temperature Annual Range (BIO5-BIO6)",
"Mean Temperature of Wettest Quarter",
"Mean Temperature of Driest Quarter",
"Mean Temperature of Warmest Quarter",
"Mean Temperature of Coldest Quarter",
"Annual Precipitation",
"Precipitation of Wettest Month",
"Precipitation of Driest Month",
"Precipitation Seasonality (Coefficient of Variation)",
"Precipitation of Wettest Quarter",
"Precipitation of Driest Quarter",
"Precipitation of Warmest Quarter",
"Precipitation of Coldest Quarter"))

r <- r[[1:19]]
names(r) <- bios$ABBRV

# Worldclim data has a scale factor of 10, so 37 is 3.7C


coords <- data.frame(x=long_grid,y=lat_grid)

points <- SpatialPoints(coords, proj4string = r@crs)

values <- extract(r, points)
values_alt <- extract(a, points)
df <- cbind.data.frame(coordinates(points),values, ALT = values_alt)
head(df)

# Choose uncorrelated variables for analysis
cormat <- cor(df)
head(round(cormat,2))

### Write a PDF of the correlation matrix
jpeg("Heatmap_AllVars.jpeg", width=9, height=8, res=500, unit="in")
  heatmap.2(abs(cormat), scale="none", col= colorRampPalette(c( "white", "blue"))(100),
          trace = "none", density.info = "none", key.xlab="|Cor(env. a,env. b)|")
dev.off()



### Scale environmental variables

df[,3:ncol(df)] <- scale(df[,3:ncol(df)])
head(df)


#check scaling worked
apply(df, 2, sd)
round(apply(df, 2, mean))

### Order 
df <- df[order(df$x),] # lowest long on left
df <- df[order(df$y, decreasing=TRUE),]

### Convert lat long to slim coordinates
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
df$slim_x <- range01(df$x)
df$slim_y <- range01(df$y)


### Subset to adaptive environments

vars <- c("x", "slim_x", "y", "slim_y", "MAT","PWarmQ", "MTDQ", "PWM", "MTWetQ", "PDM")

df_sub <- df[,which(colnames(df) %in% vars)]
head(df_sub)
cormat <- cor(df_sub)
round(cormat,1)

### Use these plots to check that SliM Maps look the same
#ggplot(df_sub) + geom_point(aes(x, y, color=MAT/10))
#ggplot(df_sub) + geom_point(aes(x, y, color=MTWetQ))

## Write scaled environment data to file ####
#plot(df_sub$slim_x, df_sub$x)
#plot(df_sub$slim_y, df_sub$y)
head(df_sub)
write.table(df_sub,"adaptive_env.txt", col.names=TRUE)

## Write uncorrelated non-adaptive envis to file
extravars <- c("x", "slim_x", "y", "slim_y", "ISO", "PSsd", "TSsd")
write.table(df[,which(colnames(df) %in% extravars)], "nuisance_env.txt", col.names=TRUE)



for (i in 3:ncol(df_sub)){
  mat <- matrix(scale(df_sub[,i]), gridsize, gridsize, byrow=TRUE)
  rownames(mat) <- rev(lats)
  colnames(mat) <- (lons)
  
  # sanity check looks good
  # head(mat[,1:5])
  # head(df_sub)
  # 
  # tail(mat[,355:360])
  # tail(df_sub)
  m2 <- t(mat)
 #don't ask me why but this produces the correct orientation in SLIM
  write.table(m2, file = paste0( colnames(df_sub)[i], "_BC_360x360.csv"), sep=",", row.names = F, col.names = F)
  
}
write.table(bios, file="bioclim.txt", col.names=TRUE)

