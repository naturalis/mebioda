
gc()
#Packages
library(sp)
library(rgdal)
library(raster)
library(biomod2)

#species used: Rhinolophus_euryale_csv2
#wordclim version: 2.0
#future variables scenario:HadGEM2ES_RCP85

#checks where your R is currently working
getwd()
#sets a new working directory
setwd("C:/Practical/")

# Load species occurrence file
sp <- read.csv("./Occurrences/corrected_data.csv",header=T) #load csv of occurrence (csv2 is Europeanstyle)
head(sp) #check table looks correct
sp_shp <- sp #rename table
coordinates(sp_shp) <- ~longitude+latitude #convert table to points shapefile
proj4string(sp_shp) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")


#Create bounding box around points
bbox <- extent(sp_shp) #create bounding box of points
bbox <- bbox+2 #increase border so we do not truncate data

#Create Binary map based on Threshold

pres <- raster("./MaxEnt/Results/Rhinolophus_euryale_Present_WLD_avg.asc")
# load present global distribution
plot(pres) #visualize

pres_clip <- crop(pres,bbox) #clip to training area
plot(pres_clip) #visualize


fut <- raster("./MaxEnt/Results/Rhinolophus_euryale_he45bi50_WLD_avg.asc")
# load future global distribution
plot(fut)
fut_clip <- crop(fut,bbox) #clip to training area
plot(fut_clip) #visualize


th <- 0.362 #define threshold
m <- c(0, th,0, th, 1, 1) #matrix everything before th as 0 and everything after th as 1
bin_mat <- matrix(m, ncol=3, byrow=TRUE) #convert to correct matrix
pres_bin <- reclassify(pres2, bin_mat) #reclassify by matrix present
fut_bin <- reclassify(fut, bin_mat) #reclassify by matrix future
plot(pres_bin)
plot(fut_bin)

#Calculate Change in Suitable Climate Conditions
range_change<-BIOMOD_RangeSize(CurrentPred=pres_bin,FutureProj=fut_bin,SpChange.Save=NULL)
#calculate range change between two maps
col.lst <- c("red3" , "gold", "grey89" , "green4") #define plot colours
plot(range_change$Diff.By.Pixel,col=col.lst) #plot range change

#Calculate occurrence area change maps
pres_bin_clip <- crop(pres_bin,bbox) #clip binary map
fut_bin_clip <- crop(fut_bin,bbox) #clip binary map
range_change_clip<-BIOMOD_RangeSize(CurrentPred=pres_bin_clip,FutureProj=fut_bin_clip,SpChange.Save=NULL)
#calculate range change between two maps
col.lst <- c("red3" , "gold", "grey89" , "green4") #define plot colours
plot(range_change_clip$Diff.By.Pixel,col=col.lst) #plot range change
plot(sp_shp,add=T, col="black", pch=19) #show points on top
