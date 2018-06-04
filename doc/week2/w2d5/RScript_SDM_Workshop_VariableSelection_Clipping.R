#############################################################################################################
#############################################################################################################
##############                                                                                 ############## 
##############                    SDM Workshop - R Script                                      ##############
##############                    Clipping + Variable Selection                                ##############
##############                    Friday 8th December 2017                                     ##############
##############                    Leon Marshall - Naturalis Biodiversity Center                ##############
##############                                                                                 ##############
#############################################################################################################
#############################################################################################################

#Packages
library(sp)
library(rgdal)
library(raster)
library(biomod2)
#functions
cor.prob <- function (X, dfr = nrow(X) - 2) {
  R <- cor(X, use="pairwise.complete.obs",  method ="spearman")
  above <- row(R) < col(R)
  r2 <- R[above]^2
  Fstat <- r2 * dfr/(1 - r2)
  R[above] <- 1 - pf(Fstat, 1, dfr)
  R[row(R) == col(R)] <- NA
  R
}
# Use this to dump the cor.prob output to a 4 column matrix
flattenSquareMatrix <- function(m) {
  if( (class(m) != "matrix") | (nrow(m) != ncol(m))) stop("Must be a square matrix.") 
  if(!identical(rownames(m), colnames(m))) stop("Row and column names must be equal.")
  ut <- upper.tri(m)
  data.frame(i = rownames(m)[row(m)[ut]],
             j = rownames(m)[col(m)[ut]],
             cor=t(m)[ut],
             p=m[ut])
}


# Load species occurrence file
sp <- read.csv("G:/SDM_Course/MAXENT/Species/Rhinolophus_euryale.csv",header=T) #load csv of occurrence
head(sp) #check table looks correct
sp_shp <- sp #rename table
coordinates(sp_shp) <- ~longitude+latitude #convert table to points shapefile
proj4string(sp_shp) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") #define projection WGS1984

#Create bounding box around points
bbox <- extent(sp_shp) #create bounding box of points
bbox <- bbox+2 #increase border so we do not truncate data
plot(bbox, col='blue') #check if box surrounds points
plot(sp_shp,add=T,pch=19,col='red') #add points

#Load Bioclim rasters
bioclims <- stack(list.files("G:/SDM_Course/Bioclims", 
                             pattern="*.bil$", full.names=TRUE, 
                                      ignore.case=TRUE)) #create a stack of present bioclims
bioclims_fut <- stack(list.files("G:/SDM_Course/Bioclims_Future", 
                             pattern="*.tif$", full.names=TRUE, 
                             ignore.case=TRUE)) #create a stack of future bioclims

#Clip rasters by bounding box
bioclims_bbox <- stack() #create empty raster stack
for(i in 1:19){
    ras <- crop(bioclims[[i]],bbox) #clip rasters by bbox
    bioclims_bbox <- stack(bioclims_bbox,ras) #stack rasters together
              } 
plot(bioclims_bbox)

#Variable Selection
bioclims_corr <- as.data.frame(bioclims_bbox, xy=T) #convert raster stack to a table
bioclims_corr <- na.omit(bioclims_corr[3:21]) #remove NA values (ocean) and remove coordiantes

sc  <- flattenSquareMatrix(cor.prob(bioclims_corr)) #Table summary correlations
sc1 <- sc[with(sc, order(cor)), ]
sc2 <- sc1[sc1$cor>=0.7|sc1$cor<=-0.7,] #High Correlations

sc2[with(sc2, order(i)), ] #list correlated variables

bioclims_corr2 <- subset(bioclims_corr, select=-c(bio18,bio14,bio17,bio3,bio5,bio10,bio6,
                                                  bio9,bio11,bio12,bio7,bio16,bio19))
#remove unwanted variables
  
sc  <- flattenSquareMatrix(cor.prob(bioclims_corr2)) #Table summary correlations
sc1 <- sc[with(sc, order(cor)), ]
sc2 <- sc1[sc1$cor>=0.7|sc1$cor<=-0.7,] #High Correlations

sc2[with(sc2, order(i)), ] #list correlated variables


#Make selection of chosen variables
list.names <- c("bio1","bio10","bio11","bio12","bio13","bio14",
                "bio15","bio16","bio17","bio18","bio19","bio2",
                "bio3","bio4","bio5","bio6","bio7","bio8","bio9") #list of names in same order as rasterstacks
select.names <- names(bioclims_corr2) #list of selected rasters
numbers <- match(select.names,list.names) #numbers corresponding to chosen variables

pres_bioclims <- bioclims[[numbers]] #select only variables interested in
fut_bioclims <- bioclims_fut[[numbers]] #select only variables interested in
clip_bioclims <- bioclims_bbox[[numbers]] #select only variables interested in

#Save bioclims as .asc
dir.create("G:/SDM_Course/MAXENT/ClimatePresent") #create folder
dir.create("G:/SDM_Course/MAXENT/ClimateFuture") #create folder
dir.create("G:/SDM_Course/MAXENT/ClimateTraining") #create folder

for(i in 1:nlayers(pres_bioclims)){
  ras_sc=scale(pres_bioclims[[1]],center=TRUE, scale=TRUE) #scale variables
  writeRaster(ras_sc,
              file=paste("G:/SDM_Course/MAXENT/ClimatePresent/",select.names[i],".asc",sep=""))
}

for(i in 1:nlayers(fut_bioclims)){
  ras_sc=scale(fut_bioclims[[1]],center=TRUE, scale=TRUE) #scale variables
  writeRaster(ras_sc,
              file=paste("G:/SDM_Course/MAXENT/ClimateFuture/",select.names[i],".asc",sep=""))
}

for(i in 1:nlayers(clip_bioclims)){
  ras_sc=scale(clip_bioclims[[1]],center=TRUE, scale=TRUE) #scale variables
  writeRaster(ras_sc, 
              file=paste("G:/SDM_Course/MAXENT/ClimateTraining/",select.names[i],".asc",sep=""))
}

#Create Binary map based on Threshold

pres <- raster("G:/SDM_Course/MAXENT/Results/Rhinolophus_euryale_ClimatePresent.asc")
# load present global distribution
plot(pres) #visualize

pres_clip <- crop(pres,bbox) #clip to training area
plot(pres_clip) #visualize


fut <- raster("G:/SDM_Course/MAXENT/Results/Rhinolophus_euryale_ClimateFuture.asc")
# load future global distribution
plot(fut)
fut_clip <- crop(fut,bbox) #clip to training area
plot(fut_clip) #visualize


th <- 0.387 #define threshold
m <- c(0, th,0, th, 1, 1) #matrix everything before th as 0 and everything after th as 1
bin_mat <- matrix(m, ncol=3, byrow=TRUE) #convert to correct matrix
pres_bin <- reclassify(pres, bin_mat) #reclassify by matrix present
fut_bin <- reclassify(fut, bin_mat) #reclassify by matrix future
plot(pres_bin)
plot(fut_bin)
#Calculate Change in Sutiable Climate Conditions
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
