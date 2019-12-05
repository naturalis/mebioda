
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

#lists all file names and paths to the files, saves it to a list
rst.names <- list.files("./Present/",pattern=".bil")[c(1,12:19,2:11)]
rst.fld <-  list.files("./Present/",pattern=".bil",full.names = T)[c(1,12:19,2:11)]
#explore your list to see if you only have rasters on it
rst.fld
#loads all rasters into a stack object (a multi-dimensional raster)
stk.present <- stack(rst.fld)

#loading future variables is a bit more complicated due to the subfolders
path.fut.26 <- list.files("./Future/he26bi50/",pattern=".tif",full.names = T)
path.fut.45 <- list.files("./Future/he45bi50/",pattern=".tif",full.names = T)
path.fut.60 <- list.files("./Future/he60bi50/",pattern=".tif",full.names = T)
path.fut.85 <- list.files("./Future/he85bi50/",pattern=".tif",full.names = T)

stk.fut.26 <- stack(path.fut.26)
stk.fut.45 <- stack(path.fut.45)
stk.fut.60 <- stack(path.fut.60)
stk.fut.85 <- stack(path.fut.85)

#checks the order of the layers loaded and renames them to bioXX, each XX represents
#a bioclimatic layer
names(stk.present)
names(stk.present) <- c("Bio01","Bio02","Bio03","Bio04",
                        "Bio05","Bio06","Bio07","Bio08",
                        "Bio09","Bio10","Bio11","Bio12",
                        "Bio13","Bio14","Bio15","Bio16",
                        "Bio17","Bio18","Bio19")

#this shows us the order is "messed up", so we need to fix it. 
names(stk.fut.26)
names(stk.fut.45)
names(stk.fut.60)
names(stk.fut.85)
path.fut.26
#the easiest way to  is just to re-load the variables again with the proper order
path.fut.26 <- list.files("./Future/he26bi50/",pattern=".tif",full.names = T)[c(1,12:19,2:11)]
path.fut.45 <- list.files("./Future/he45bi50/",pattern=".tif",full.names = T)[c(1,12:19,2:11)]
path.fut.60 <- list.files("./Future/he60bi50/",pattern=".tif",full.names = T)[c(1,12:19,2:11)]
path.fut.85 <- list.files("./Future/he85bi50/",pattern=".tif",full.names = T)[c(1,12:19,2:11)]

stk.fut.26 <- stack(path.fut.26)
stk.fut.45 <- stack(path.fut.45)
stk.fut.60 <- stack(path.fut.60)
stk.fut.85 <- stack(path.fut.85)

#check if they loaded fine
names(stk.fut.26)
names(stk.fut.45)
names(stk.fut.60)
names(stk.fut.85)

#and then rename them easily
list.of.names <- c("Bio01","Bio02","Bio03","Bio04",
                   "Bio05","Bio06","Bio07","Bio08",
                   "Bio09","Bio10","Bio11","Bio12",
                   "Bio13","Bio14","Bio15","Bio16",
                   "Bio17","Bio18","Bio19")

names(stk.fut.26) <- list.of.names
names(stk.fut.45) <- list.of.names
names(stk.fut.60) <- list.of.names
names(stk.fut.85) <- list.of.names

names(stk.fut.26)
# Load species occurrence file
# notice im using read.csv2, which expects a EU type of table. If you want to use the NA style,
#then you must switch read.csv2 with read.csv
#You can also use custom delimitrs 
sp <- read.csv2("./Occurrences/Rhinolophus_euryale_csv2.csv",header=T) #load csv of occurrence
head(sp) #check table looks correct
sp_shp <- sp #rename table
coordinates(sp_shp) <- ~longitude+latitude #convert table to points shapefile
proj4string(sp_shp) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")


#Create bounding box around points
bbox <- extent(sp_shp) #create bounding box of points
bbox <- bbox+2 #increase border so we do not truncate data
plot(stk.present$Bio01,ext=bbox+2)
plot(bbox, col='blue',add=T) #check if box surrounds points
plot(sp_shp,add=T,pch=19,col='red') #add points

#cropping the present data
stk.present.AOI.crop <- crop(stk.present,bbox) #clip to training area
#plotting the example
par(mfrow=c(1,2)) #sets the plotting area to a 1 line 2 columns set up
plot(stk.present$Bio01,main="Original extent")
plot(stk.present.AOI.crop $Bio01,main="Cropped extent")
par(mfrow=c(1,1)) #sets it back to 1 image per plot

#cropping the future data
stk.fut.26.AOI.crop <- crop(stk.fut.26,bbox)
stk.fut.45.AOI.crop <- crop(stk.fut.45,bbox)
stk.fut.60.AOI.crop <- crop(stk.fut.60,bbox)
stk.fut.85.AOI.crop <- crop(stk.fut.85,bbox)

#you can see the 
par(mfrow=c(1,2))
plot(stk.present$Bio01, main="Present data")
plot(stk.fut.26$Bio01,main="Future data")


#cropping the present world data
stk.present <- crop(stk.present,stk.fut.26)

#now we can save them to another folder in a format
#that maxent can read
#saving the cropped present data in .asc format
writeRaster(stk.present.AOI.crop,
            "./Present_AOI/.asc",
            overwrite=T,
            bylayer=T,
            suffix="names")

#saving present data uncropped - notice, these files might be extra large
writeRaster(stk.present,
            "./Present_WLD/.asc",
            overwrite=T,
            bylayer=T,
            suffix="names")

#saving scenario he26bi50 uncropped and cropped
writeRaster(stk.fut.26,
            "./Future_WLD/he26bi50_WLD/.asc",
            overwrite=T,
            bylayer=T,
            suffix="names")

writeRaster(stk.fut.26.AOI.crop,
            "./Future_AOI/he26bi50_AOI/.asc",
            overwrite=T,
            bylayer=T,
            suffix="names")

#saving scenario he45bi50 uncropped and cropped
writeRaster(stk.fut.45,
            "./Future_WLD/he45bi50_WLD/.asc",
            overwrite=T,
            bylayer=T,
            suffix="names")

writeRaster(stk.fut.45.AOI.crop,
            "./Future_AOI/he45bi50_AOI/.asc",
            overwrite=T,
            bylayer=T,
            suffix="names")

#saving scenario he50bi50 uncropped and cropped
writeRaster(stk.fut.60,
            "./Future_WLD/he60bi50_WLD/.asc",
            overwrite=T,
            bylayer=T,
            suffix="names")

writeRaster(stk.fut.60.AOI.crop,
            "./Future_AOI/he60bi50_AOI/.asc",
            overwrite=T,
            bylayer=T,
            suffix="names")

#saving scenario he85bi50 uncropped and cropped
writeRaster(stk.fut.85,
            "./Future_WLD/he85bi50_WLD/.asc",
            overwrite=T,
            bylayer=T,
            suffix="names")

writeRaster(stk.fut.85.AOI.crop,
            "./Future_AOI/he85bi50_AOI/.asc",
            overwrite=T,
            bylayer=T,
            suffix="names")

#once the last step is done, you can proceed to the next part to 
#test for autocorrelation between variables.
#consider that if you close R meanwhile, you have to reload all the variables
#again

