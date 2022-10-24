
# Spatial area

# SW 49.923502, -123.191617
# NE 53.407794, -114.808504

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
names(r) <- bios$DESC

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
max(abs(lower.tri(cormat)))
heatmap(abs(cormat), scale="none")

vars <- c("x", "y", "MAT","PWarmQ", "MTDQ", "PWM", "MTWetQ", "PDM")

df_sub <- df[,which(colnames(df) %in% vars)]
head(df_sub)
cormat <- cor(df_sub)
round(cormat,1)
heatmap(abs(cormat), scale="none")

df_sub <- df_sub[order(df_sub$x),] # lowest long on left
df_sub <- df_sub[order(df_sub$y, decreasing=TRUE),]
head(df_sub)

ggplot(df_sub) + geom_point(aes(x, y, color=MAT/10))
ggplot(df_sub) + geom_point(aes(x, y, color=MTWetQ))

df_sub[,3:ncol(df_sub)] <- scale(df_sub[,3:ncol(df_sub)])

#check scaling worked
apply(df_sub, 2, sd)
apply(df_sub, 2, mean)

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
  write.table(m2, file = paste0("bioclim/", colnames(df_sub)[i], "_BC_360x360.csv"), sep=",", row.names = F, col.names = F)
  
}
write.table(bios, file="bioclim/bioclim.txt", col.names=TRUE)

