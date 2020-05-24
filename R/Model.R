##########################
#     Load libraries     #
##########################

library(spocc)
library(spThin)
library(dismo)
library(rgeos)
library(ENMeval)
library(dplyr)
library(rgdal)
library(rJava)
library(maptools)
library(ggplot2)
library(devtools)
library(ENMGadgets)
library(corrplot)
library(ntbox)
library(wallace)

########################
#      Load data       #
########################

setwd("~/Documentos/Omar/Projects/RangeExtentions/E.helias/")

source("../R/thresh.R")

##########################################
#LoadPoints
occ <- read.csv("All_occ_2.csv", header = T)

#Quickly Occ Clean
occs.dups <- duplicated(occ[c('Lon', 'Lat')])
occ <- occ[!occs.dups,]
occ
occ$occID <- row.names(occ)

#########################################
#Load M polygon
Mpol <- readOGR("shp/M3_base.shp")

#########################################
#Read WorldClim variables
bio1 <- raster("~/Documentos/Omar/Qgis/WorldClim/Bio_2.5m/wc2.1_2.5m_bio_1.tif")
bio3 <- raster("~/Documentos/Omar/Qgis/WorldClim/Bio_2.5m/wc2.1_2.5m_bio_3.tif")
bio4 <- raster("~/Documentos/Omar/Qgis/WorldClim/Bio_2.5m/wc2.1_2.5m_bio_4.tif")
bio5 <- raster("~/Documentos/Omar/Qgis/WorldClim/Bio_2.5m/wc2.1_2.5m_bio_5.tif")
bio6 <- raster("~/Documentos/Omar/Qgis/WorldClim/Bio_2.5m/wc2.1_2.5m_bio_6.tif")
bio7 <- raster("~/Documentos/Omar/Qgis/WorldClim/Bio_2.5m/wc2.1_2.5m_bio_7.tif")
bio10 <- raster("~/Documentos/Omar/Qgis/WorldClim/Bio_2.5m/wc2.1_2.5m_bio_10.tif")
bio11 <- raster("~/Documentos/Omar/Qgis/WorldClim/Bio_2.5m/wc2.1_2.5m_bio_11.tif")
bio12 <- raster("~/Documentos/Omar/Qgis/WorldClim/Bio_2.5m/wc2.1_2.5m_bio_12.tif")
bio13 <- raster("~/Documentos/Omar/Qgis/WorldClim/Bio_2.5m/wc2.1_2.5m_bio_13.tif")
bio16 <- raster("~/Documentos/Omar/Qgis/WorldClim/Bio_2.5m/wc2.1_2.5m_bio_16.tif")

envs <- stack(bio1,bio3,bio4,bio5,bio6,bio7,
              bio10,bio11,bio12,bio13,bio16)


###########################################
#            Reduce Sampling Bias         #
###########################################

occ2 <- occ
coordinates(occ2) <- ~Lon+Lat ; proj4string(occ2) <- CRS(proj4string(Mpol))

ras <- raster(occ2) #Create raster layer
res(ras) <- 1 #Set grid resolution 1°
ras <- extend(ras, extent(ras)+1) #Expand raster layer

occS <- gridSample(occ2, ras, n = 1) #Just one occ per cell

occ2 <- as.data.frame(coordinates(occS))
colnames(occ2) <- c("Lon", "Lat")

occ2 <- occ2
coordinates(occ2) <- ~Lon+Lat
proj4string(occ2) <- CRS(proj4string(Mpol))

###########################################
#      Get records whitin the M polygon   #
###########################################

inOcc <- over(occ2, Mpol)

occ2 <- as.data.frame(coordinates(occ2))
colnames(occ2) <- c("Lon","Lat")

occ3 <- occ2[-which(is.na(inOcc)), ]
occ3 <- as.data.frame(coordinates(occ3))
colnames(occ3) <- c("Lon","Lat")

occ4 <- occ2[which(is.na(inOcc)),]
occ4 <- as.data.frame(coordinates(occ4))
colnames(occ4) <- c("Lon","Lat")

plot(Mpol)
points(occ2$Lon, occ2$Lat, col="red")
points(occ3$Lon, occ3$Lat, col="green")
points(occ4$Lon, occ4$Lat, col="blue")

#################################################
#   Remove occ without environmental values     #
#################################################

# extract environmental values for each occurrence
locs.vals <- raster::extract(envs[[1]], occ3[, c('Lon', 'Lat')])
# remove occs without environmental values
occ3 <- occ3[!is.na(locs.vals), ] 

#################################################
#           Create Background data              #
#################################################

#Buffer size of the study extent polygon defined as 1 degrees.
bgExt <- Mpol
bgExt <- rgeos::gBuffer(bgExt, width = 1)

# crop the environmental rasters by the background extent shape
  envsBgCrop <- raster::crop(envs, bgExt)
envsBgMsk <- raster::mask(envsBgCrop, bgExt)

#save(envsBgMsk, file="envsBgMsk.rda")
#load("envsBgMsk.rda")

# sample random background points
bg.xy <- dismo::randomPoints(envsBgMsk, 20000)
# convert matrix output to data frame
bg.xy <- as.data.frame(bg.xy)  

plot(bgExt)
points(bg.xy)
points(occ3$Lon, occ3$Lat, col="red")

##########################################################
#                Create partitioning Data                #
##########################################################

#Occ3 was preoviosly created !!
occs.xy <- occ3[, c('Lon','Lat')]

#Create group partition
group.data <- ENMeval::get.block(occ=occs.xy, bg.coords=bg.xy)

occs.grp <- group.data[[1]]
bg.grp <- group.data[[2]]

par(mfrow=c(1,2))
plot(bgExt, main="Original Occ")
points(occs.xy, col=occs.grp)
plot(bgExt, main="Background Occ")
points(bg.xy, col=bg.grp)
par(mfrow=c(1,1))

##########################################################
#                 Build and Evaluate ENM                 #
##########################################################

# define the vector of regularization multipliers to test
rms <- seq(1, 5, .5)

# iterate model building over all chosen parameter settings
e <- ENMeval::ENMevaluate(occs.xy, envsBgMsk, bg.coords = bg.xy, RMvalues = rms, fc = c('L', 'LQ', 'LQH', 'LQHP'), 
                          method = 'user', occs.grp, bg.grp, clamp = TRUE, algorithm = "maxent.jar")

#save(e, file = "Maxent/Models/LQHP_3/e.rda")

##########################################################
#            Extract results from the analysis           #
##########################################################

evalMods <- e@models
names(evalMods) <- e@results$settings
evalTbl <- e@results
evalPreds <- e@predictions

##########################################################
#                 Load the G area polygon                #
##########################################################

projPoly <- readOGR("shp/Proj_pol.shp")
predsProj <- raster::crop(envs, projPoly)
predsProj <- raster::mask(predsProj, projPoly)

#save(predsProj, file="predsProj.rda")
#load("predsProj.rda")

##########################################################
#                  See model plots                       #
##########################################################

ENMeval::eval.plot(evalTbl, value = "avg.test.orMTP")
ENMeval::eval.plot(evalTbl, value = "avg.test.or10pct")
ENMeval::eval.plot(evalTbl, value = "AICc")
ENMeval::eval.plot(evalTbl, value = "avg.test.AUC")

evalTbl

#write.csv(evalTbl, "Maxent/evalTbl_2.5m.csv", quote = F, row.names = F)

#str(evalMods$LQ_1)

#evalMods$H_3@html

##########################################################
#               Model Selection Rules                    #
##########################################################

MTPtable <- evalTbl[which(evalTbl$avg.test.orMTP<0.05),]
MTPtable

AIC <- MTPtable$settings[which(MTPtable$delta.AICc == min(MTPtable$delta.AICc))]
AIC

AUC <- MTPtable$settings[which(MTPtable$avg.test.AUC == max(MTPtable$avg.test.AUC))]
AUC

#AIC LQHP_5
#AUC LQ_3

evalMods$LQHP_3@path

##########################################################
#                  Predict over G                        #
##########################################################

pred <- dismo::predict(evalMods$LQHP_3, envsBgMsk, args="outputformat=cloglog")
plot(pred)

#proc <- pROC(continuous_mod = pred, test_data = test.data, n_iter = 100)
#proc

pred <- dismo::predict(evalMods$LQHP_3, predsProj, args="outputformat=cloglog")
plot(pred)

#proc <- pROC(continuous_mod = pred, test_data = occ4, n_iter = 100)
#proc

#########################################################
#                  Threshold and binary                 #
#########################################################

occPredVals <- raster::extract(pred, occs.xy)

thr <- thresh(modOccVals = occPredVals, type = "p10")
thr

pred.bin <- pred>thr
plot(pred.bin)

writeRaster(pred, filename = "Maxent/Models/LQHP_3/Pred_LQHP3_AIC", format = "GTiff")
writeRaster(pred.bin, filename = "Maxent/Models/LQHP_3/Pred_LQHP3_AICbin", format ="GTiff")
