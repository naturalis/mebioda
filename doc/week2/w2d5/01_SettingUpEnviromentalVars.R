
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
rst.names <- list.files("./Present/",pattern=".tif")
rst.fld <-  list.files("./Present/",pattern=".tif",full.names = T)
#explore your list to see if you only have rasters on it
rst.fld
#loads all rasters into a stack object (a multi-dimensional raster)
my.stack <- stack(rst.fld)

names(my.stack)
names(my.stack) <- c("Bio01","Bio02","Bio03","Bio04",
                     "Bio05","Bio06","Bio07","Bio08",
                     "Bio09","Bio10","Bio11","Bio12",
                     "Bio13","Bio14","Bio15","Bio16",
                     "Bio17","Bio18","Bio19")

# Load species occurrence file
sp <- read.csv2("./Occurrences/Rhinolophus_euryale_csv2.csv",header=T) #load csv of occurrence
head(sp) #check table looks correct
sp_shp <- sp #rename table
coordinates(sp_shp) <- ~longitude+latitude #convert table to points shapefile
proj4string(sp_shp) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")


#Create bounding box around points
bbox <- extent(sp_shp) #create bounding box of points
bbox <- bbox+2 #increase border so we do not truncate data
plot(my.stack$Bio01,ext=bbox+2)
plot(bbox, col='blue',add=T) #check if box surrounds points
plot(sp_shp,add=T,pch=19,col='red') #add points


#cropping the rasters
rst.AOI.crop <- crop(my.stack,bbox) #clip to training area
par(mfrow=c(1,2))
plot(my.stack$Bio01,main="Original extent")
plot(rst.AOI.crop$Bio01,main="Cropped extent")

#now we can save them to another folder in a format
#that maxent can read
writeRaster(rst.AOI.crop,
            "./Present_rdy/.asc",
            overwrite=T,
            bylayer=T,
            suffix="names")


### pairwise testing 

#first we convert the cropped raster to a data.frame
df.crop.stack <- na.omit(as.data.frame(rst.AOI.crop)) #we also remove NA's
#now this stores the pearson correlation in a matrix
cor.tab <-cor(df.crop.stack)
write.csv2(cor.tab,"PearsonR_Europe.csv")


#multicollinearity testing
library(usdm)

#e.g. i select Bio01; Bio04; Bio07; Bio 12; Bio 15 and bio 19
head(df.crop.stack)
df.crop.stack.selection <- df.crop.stack[,c(1,4,7,12,15,19)] #select only the variables i am interested on
head(df.crop.stack.selection)

vif(df.crop.stack.selection, maxobservations=nrow(df.crop.stack.selection))

df.crop.stack.selection <- df.crop.stack[,c(1,4,12,15,19)] #minos the temperature range
vif(df.crop.stack.selection, maxobservations=nrow(df.crop.stack.selection))



#Last step: 

#Saving the cropped variables to the correct folder, in the .asc format

#we had already saved the cropped rast
path2rst.HadGEM2ES_RCP85 <- list.files("./Future/", full.names = T)
path2rst.HadGEM2ES_RCP85
#in this case the order of variables is changed, its better to keep everything in the
#same order. It's easy to adapt the code to read in our preferred order
#WARNING: this step might be different in your case, CHECK IT FIRST
rst.fut.stack <- stack(path2rst.HadGEM2ES_RCP85[c(1,12:19,2:11)])
rst.fut.stack.crop <- crop(rst.fut.stack,bbox)
names(rst.fut.stack.crop)

#now we can jus save the variables to a new folder
#but first, we rename the variables to the same as the present vars
names(rst.fut.stack.crop)<- names(rst.AOI.crop)

writeRaster(rst.fut.stack.crop,
            "./Future_rdy/.asc",
            overwrite=T,
            bylayer=T,
            suffix="names")

#now you can go to maxent




#correcting occurrences
sp <- read.csv2("./Occurrences/Rhinolophus_euryale_csv2.csv",header=T) #load csv of occurrence
unique(sp$species)
write.csv(sp,"./Occurrences/Rhinolophus_euryale_csv1.csv",row.names = F)
