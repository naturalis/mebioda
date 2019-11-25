#Installing packages - online resources

#Firt steps
      #route one - Using the console
      install.packages ("raster") #Copy this to console and press enter
      
      #if all goes well, the console will say: (some extra packages will be installed because it is the first time raster is installed)
        #package ‘sp’ successfully unpacked and MD5 sums checked
        #package ‘Rcpp’ successfully unpacked and MD5 sums checked
        #package ‘raster’ successfully unpacked and MD5 sums checked
      
      #route one - Using the script 
      install.packages ("raster") #highlight this section and press ctrl + enter
      
        #package ‘raster’ successfully unpacked and MD5 sums checked
      
      #using the packages section of the user interface:
      
        #Press Install
        #Check the repository:
          #CRAN - online repository, can be changed in other options to access private repositories
          #Package - You can download the package and point R studio to install it
      
        #write package name (R Studio will auto-complete)
        #Press enter...

      
#Second steps
      
        #These custom options available through the R studio interface are obviously available in the previous console commands. 
          
          #https://www.rdocumentation.org/packages/utils/versions/3.5.1/topics/install.packages
          #using the console command allows you to install multiple packages with one command
      
        #Run this!
        install.packages(c("rgdal","sp","maptools"))
        
          #What is the c for?? 
          # https://stackoverflow.com/questions/11488820/why-use-c-to-define-vector
          # R stores numbers as vectors, arrays, lists, etc (typical object types in programming)
          # my.group <- c(vecto1, vector2, ...., vectorn) 
          # my.group is now a vector of type (etc) that holds all the data

        #Run this
        c(1,3,4,5,"ola")
        typeof(c(1,3,4,5,"ola"))
        c(1,3,4,5)
        typeof(c(1,3,4,5))
        
        #When you provide the input c("rgdal","sp","maptools") to install.packages function you are telling R to first install package 1, then iterate through
        
#installing useful packages for this course.
        
        #here is a list of very useful packages - install them the way you prefer
        raster
        sp
        ggmap
        rgdal
        rgeos
        maptools
        dplyr + tidyr
        dismo
        biomod2 
		RStoolbox

        #extra stuff that may or may not be useful
        
        #accuraccy and models evaluation
        PresenceAbsence
        psych
        e1071
        caret
        
#How to load an installed package?
        
        library(raster)
        
        #> library(raster)
        #Loading required package: sp
        #>
        
        #what if the package is not installed?
        library(dismo)
        #> library(dismo)
        #Error in library(dismo) : there is no package called ‘dismo’
        #> 
        
#How to unload a package? .. but why would you do that? :P
        
        detach("package:raster", unload=TRUE)
         
        
#what happens when a package is not installed/loaded and you try to run a function?
        
        
        #downloading worldclim data -------------------------
        w.stack = getData('worldclim',
                          var='bio',
                          res=10) #Downloading a WCLIM tile set

        #Error in getData("worldclim", var = "bio", res = 10) : could not find function "getData"
        
        #Sometimes it even tells which package is missing. How to solve? load the package. If you don't know which package it is.. google
        library(raster)
        w.stack = getData('worldclim',
                          var='bio',
                          res=10) #Downloading a WCLIM tile set