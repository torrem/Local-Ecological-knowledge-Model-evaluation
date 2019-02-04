# Produce SI plots for variables in HSI Model
# HSI Stat Areas 467,511,512,513,514,515

library(classInt)
library(mgcv)
library(fmsb)
library(maptools)
library(rgdal)
library(gstat)
library(automap)
library(visreg)
library(colorRamps)
library(plyr)
library(geosphere)
library(raster)
library(readr)
library(rgeos) 
library(spatstat)
library(dplyr)
library(randomForest)


options(scipen = 999)


# word.tif = function(filename="Word_Figure_%03d.tif", zoom=3, width=20, height=20, pointsize=14, ...) {
#   if (!grepl("[.]ti[f]+$", filename, ignore.case=TRUE))
#     filename = paste0(filename,".tif")
#   tiff(filename=filename, compression="lzw", res=600*zoom,
#        width=width, height=height, units='cm', pointsize=pointsize, ...)
# }
# 
# #font to Times
# windowsFonts(Times=windowsFont("TT Times New Roman"))



setwd("C:/Users/Mike/Documents/Scallop/Data/LEK")

COAST = readShapePoly("C:/Users/Mike/Documents/Scallop/Data/GIS/US_Coastlines/US_Coast.shp")




scallop = read.csv("C:/Users/Mike/Documents/Scallop/Data/LEK/Cob ScallopALT.csv")
scallop = subset(scallop, Sediment18 =="G"|Sediment18 =="S"|Sediment18 =="M"|Sediment18 =="R")
scallop$surveydept = as.numeric(scallop$surveydept)

#Cob = read.csv("C:/Users/Mike/Documents/Scallop/Data/LEK/CobPoints.csv")
Cob = read.csv("C:/Users/Mike/Documents/Scallop/Data/LEK/CobPoints2.csv")

Cob$depth = Cob$depth*-1
#Cob = subset(Cob, depth > 0)
#plot(Cob$Lat ~ Cob$Long)


scallop = scallop[,c(5,6,9,14,19)]
scallop = droplevels(scallop)
scallop = subset(scallop, abundance<2000)



split=sample(1:(nrow(scallop)), size = nrow(scallop)*0.90)
test = scallop[split,]
train = scallop[-split,]


#Cob = Cob[,2:5]
names(Cob)[1]="long_start"
names(Cob)[2]="lat_start"
names(Cob)[3]="Sediment18"
names(Cob)[4] = "surveydept"




 rf=randomForest(abundance ~    surveydept + Sediment18+ long_start + lat_start ,
                 mtry=1.1,data = train)
 rf


 
 gm = gam(abundance~s(surveydept, bs="ts") + s(lat_start,bs="ts")+ Sediment18  
         + s(long_start,bs="ts")+s(lat_start,long_start,bs="ts")
         + s(lat_start,surveydept,bs="ts") +s(surveydept,long_start,bs="ts"), 
         data=train, family = "tw(theta = NULL, link = 'log', a=1.01, b=1.99)", select=T )
 summary(gm)
 
 
 

 
 Cob = cbind(Cob, predabu = predict(rf, newdata = Cob, type='response'))
 
 e = extent(min(Cob$long_start), max(Cob$long_start), min(Cob$lat_start), max(Cob$lat_start)    )
 tr=105
 nr = tr*1.01
 nc = as.numeric(tr*(sqrt((e[1] - e[2])^2)/sqrt((e[3] - e[4])^2)))
 r <- raster(e, ncol=nc, nrow=nr)
 
 
 Cobras <- rasterize(Cob[, c(1,2)], r, Cob[,5], fun=mean)
 # rasabuCom <- rasterize(rasdatCom[, c(1,2)], r, rasdatCom[,9], fun=mean)
 
 #par(mfrow=c(1,2))
 mypal <- colorRampPalette(c("blue", "cyan", "yellow", "red"), bias = 1)
# word.tif('Abu_plot96')
 plot(Cobras, col = mypal(25))
# dev.off()
 
 
 #writeRaster(Cobras, "Cobras", format = "ascii")
 
 
 #----------------------------Compare Model results to LEK-------------------
 
 LEK = readOGR(dsn = "C:/Users/Mike/Documents/Scallop/Data/LEK/Shapes", layer = "COBLEK")
 
 LEK <- LEK[LEK@data$Intvw_scor!=0,]
 
 coordinates(Cob)=~long_start+lat_start
 
 projection(LEK, asText=TRUE)
 projection(Cob) <- projection(LEK, asText=TRUE)
 projection(Cob, asText=TRUE)
 
 join = over(LEK, Cob[,"predabu"], FN= mean)
 
 f = cbind(data.frame(LEK),join)
 
 f = f[!is.na(f$predabu),]
 
 f$Zbin = f$Intvw_scor
 
 f$Zbin = ifelse(f$Zbin > 0, "H", "L")
 
 f$predbin = f$predabu
 

CI =  classIntervals(f$predabu, n=2, style="jenks")
 
median(f$predabu)
#mean(f$predabu)
f$predbin = ifelse(f$predbin > 176, "H", "L")
 #f$predbin = ifelse(f$predbin > CI$brks[[2]], "H", "L")
 
# f$sumER = ifelse(f$predbin == f$Zbin, 0, 1)
 
# 1-(sum(f$sumER)/nrow(f))
 
 
 
 
 
 modelPred = as.character(f$predbin);LEKScore = f$Zbin
 
 MK = data.frame(cbind(modelPred,LEKScore))
 
 a = length(which(MK$modelPred == "L" & MK$LEKScore == "L"))
 b = length(which(MK$modelPred == "H" & MK$LEKScore == "L"))
 c = length(which(MK$modelPred == "L" & MK$LEKScore == "H"))
 d = length(which(MK$modelPred == "H" & MK$LEKScore == "H"))
 
 Sensitivity = d/(c+d);Sensitivity
 Specificity = a/(a+b);Specificity
 Accuracy = (a+d)/(a+b+c+d); Accuracy
 
 
 








#------------------------END-------------------------------








# library(dismo)
# library(gbm)
# 
# 
# word.tif = function(filename="Word_Figure_%03d.tif", zoom=3, width=15, height=15, pointsize=12, ...) {
#   if (!grepl("[.]ti[f]+$", filename, ignore.case=TRUE))
#     filename = paste0(filename,".tif")
#   tiff(filename=filename, compression="lzw", res=600*zoom,
#        width=width, height=height, units='cm', pointsize=pointsize, ...)
# }
# 
# #font to Times
# windowsFonts(Times=windowsFont("TT Times New Roman"))
# 
# 
# 
# setwd("C:/Users/Mike/Documents/Scallop/Data/GAM")
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# scallop = read.csv("scallop data17DIST.csv")
# 
# scallop = subset(scallop, wq < 0.04)
# 
# 
# d=scallop
# 
# d = d[,c(11,7:8,13:16)]
# 
# d$abundance = (log(d$abundance+1))
# 
# #Identify optimal number of trees (nt)
# #start with learning rate (lr) of 0.01
# 
# lr01 = gbm.step(data = d, gbm.x = c(1,2,4,5), gbm.y = 3, 
#                 family = "gaussian", tree.complexity = 5, learning.rate = 0.04,
#                 bag.fraction = 0.5)
# 
# #look at variable importance
# 
# #word.tif('BRT')
# summary(lr01)
# #dev.off()
# 
# #interrofate and plot interactions
# #assess the extent to which pairwise interactions exist in the data
# int = gbm.interactions(lr01)
# f =int$interactions
# int$rank.list
