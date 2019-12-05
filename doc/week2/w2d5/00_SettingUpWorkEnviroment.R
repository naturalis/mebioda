
#This sets the work directory on your C drive and builds from there
setwd("C:/")
dir.create("./Practical/")
setwd("./Practical")
#from now on we are operating on C:/Practical/

#creating each folder
dir.create("downloads")
dir.create("Occurrences")

dir.create("Future")
dir.create("Future_AOI")
dir.create("Future_WLD")

dir.create("Present")
dir.create("Present_AOI")
dir.create("Present_WLD")

dir.create("Maxent")
dir.create("./Maxent/Results")
dir.create("R_Scripts")

#assuming you use the same scenario as me - Careful here
setwd("c:/Practical/Future/")
scenario.list <- c("./he26bi50/",
                   "./he45bi50/",
                   "./he60bi50/",
                   "./he85bi50/")
for (i in scenario.list){dir.create(i)}
setwd("c:/Practical/Future_AOI/")
scenario.list_AOI <- c("./he26bi50_AOI/",
                       "./he45bi50_AOI/",
                       "./he60bi50_AOI/",
                       "./he85bi50_AOI/")
for (i in scenario.list_AOI){dir.create(i)}
setwd("c:/Practical/Future_WLD/")
scenario.list_WLD <- c("./he26bi50_WLD/",
                       "./he45bi50_WLD/",
                       "./he60bi50_WLD/",
                       "./he85bi50_WLD/")
for (i in scenario.list_WLD){dir.create(i)}

