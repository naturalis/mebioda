

gc()
#Packages
library(sp)
library(rgdal)
library(raster)


###first load all the env data of your AOI in the present - you can skip if you have not closed R
### or cleared the data in R
setwd("c:/Practical")
path2presentData <- list.files("./Present_AOI/",pattern=".asc",full.names = T)
stk.present.AOI.crop <- stack(path2presentData)
names(stk.present.AOI.crop) <- c("Bio01","Bio02","Bio03","Bio04",
                                 "Bio05","Bio06","Bio07","Bio08",
                                 "Bio09","Bio10","Bio11","Bio12",
                                 "Bio13","Bio14","Bio15","Bio16",
                                 "Bio17","Bio18","Bio19")

### the autocorrelation testing is important ONLY for the areas where the
### model is trained, so, for this section, we use only the cropped enviromental data

### pairwise testing 
#first we convert the cropped raster to a data.frame
stk.present.AOI.crop <- na.omit(as.data.frame(stk.present.AOI.crop)) #we also remove NA's
#now this stores the pearson correlation in a matrix
cor.tab <-cor(stk.present.AOI.crop)
#remember to change to write.csv if needed
write.csv2(cor.tab,"CorrelationTable_AOI.csv")


#multicollinearity testing
library(usdm)

#e.g. i select Bio01; Bio04; Bio07; Bio 12; Bio 15 and bio 19
head(stk.present.AOI.crop)
#select only the variables i am interested on
df.stk.AOI <- stk.present.AOI.crop[,c(1,4,7,12,15,19)] 
#confirm the selection
head(df.stk.AOI)
#and the VIF test
vif(df.stk.AOI, maxobservations=nrow(df.stk.AOI))

df.stk.AOI <- stk.present.AOI.crop[,c(1,4,12,15,19)] #minus the temperature range
vif(df.stk.AOI, maxobservations=nrow(df.stk.AOI))


