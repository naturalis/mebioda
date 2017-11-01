#####needed libraries
library(ape)
library(apTreeshape)
library(phytools)

#setwd("C:/CUSTOM_R_FILES")  #location of input and output files
setwd("C:/R_MESADATA3")

######GLOBALS#####
SEED_VECTOR = c(4991) #,4915,4914,4913,4912,4911,4910)
SVSIZE = as.numeric(length(SEED_VECTOR))

TIME_VECTOR = c(0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150,155,160,165,170,175,180,185,195,200,205,210,215,220,225,230,235,240,245,250,255,260,265,270,275,280,285,290,295,300,301,305,310,315,320,325,330,335,340,345,350,355,360,365,370,375,380,385,390,395,400,405,410,415,420,425,430,435,440,445,450,455,460,465,470,475,480,485,490,495,500,505,510,515,520,525,530,535,540,545,550,555,560,565,570,575,580,585,590,595,600)
#TIME_VECTOR = c(0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150,155,160,165,170,175,180,185,190,195,200,205,210,215,220,225,230,235,240,245,250,255,260,265,270,275,280,285,290,295,300,301,305,310,315,320,325,330,335,340,345,350,355,360,365,370,375,380,385,390,395,400,405,410,415,420,425,430,435,440,445,450,455,460,465,470,475,480,485,490,495,500,505,510,515,520,525,530,535,540,545,550,555,560,565,570,575,580,585,590,595,600)
#TIME_VECTOR = c(0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150,155,160,165,170,175,180,185,190,195,200,205,210,215,220,225,230,235,240,245,250,255,260,265,270,275,280,285,290,295,300,301,305,310,315,320,325,330,335,340,345,350,355,360,365,370,375,380,385,390,395,400,405,410,415,420,425,430,435,440,445,450,455,460,465,470,475,480,485,490,495,500,505,510,515,520,525,530,535,540,545,550,555,560,565,570,575,580,585,590,595,600,605,610,615,620,625,630,635,640,645,650,655,660,665,670,675,680,685,690,695,700,705,710,715,720,725,730,735,740,745,750,755,760,765,770,775,780,785,790,795,800,805,810,815,820,825,830,835,840,845,850,855,860,865,870,875,880,885,890,895,900,905,910,915,920,925,930,935,940,945,950,955,960,965,970,975,980,985,990,995,1000)
TVSIZE = as.numeric(length(TIME_VECTOR))

NEXFILEEXT = '_1.nex'

SERIESSIZE = 1  #this is the number of treatment types

EXT_STRENGTH_STR = '0.9'  #the strength of extinction, e.g. 0.9 = 90% extinction

#vector containing names of data we will have in output
STATSVEC = c("TIME", "extant size", "COLL_ROUGH",
             "COLL_ROUGH_SIZESTD", "COLL_HEARD", "COLL_BLUM", "SACK_BLUM",
             "ALDOUS_B", "B_LCI", "B_UCI", "PHG", "ROOT_TIP_DIST")
SVLEN = as.numeric(length(STATSVEC))

###############end GLOBALS

######a function to find the Newick tree structure in a Nexus file and import it
GetTreeFromNexusFile = function(ITF)
{
  INFILE = file(ITF,"r")
  INLINE = readLines(INFILE,n=1)  #get 1st line to make sure not stuck
  message('Finding tree structure...')
  while(TRUE)
  {
    INLINE = readLines(INFILE,n=1)
    if (length(INLINE) == 0) {
      message('reached end of file')
      break
    }
    else {
      #message(INLINE)
      
      TEST = grep("TREE * ", INLINE)
      if (length(TEST) > 0) {
        GRABLINE = INLINE
        message('Got tree structure!')
      }
    }
  }
  close(INFILE)
  
  GRABSPLIT = strsplit(GRABLINE,split="] ")
  TREESTRUCT = GRABSPLIT[[1]][2]
  INTREE = read.tree(text=TREESTRUCT)
  #-------------------------------------------------------------#
  
  summary(INTREE)
  
  return(INTREE)
}
#####

######a function to remove all extinct leaves from a tree and return only the extant ones
DropExtinctLeaves = function(ITREE)
{
  message('Getting extincts...')
  tips = getExtinct(ITREE, tol=0.005)
  Ntip(ITREE) - length(tips)
  message('Dropping extincts...')
  ITREE.EXTANT = drop.tip(ITREE,tips)
  summary(ITREE.EXTANT)
  return (ITREE.EXTANT)
}
#####

#####a function to call the calculation of the balance statistics using apTreeshape routines
CalcBalanceStats = function(CURR_TIME, ITREE.EXTANT)
{
  RESULTS.TEMP = matrix(data=0, nrow=1, ncol=(SVLEN))
  RESULTS.TEMP[1] = CURR_TIME
  RESULTS.TEMP[2] = Ntip(ITREE.EXTANT)
  
  if (Ntip(INTREE.EXTANT) >= 3)
  {
    colnames(RESULTS.TEMP) = STATSVEC
    message('Calculating tree shape stats...')
    INTREE.EXTANT.TS = as.treeshape(INTREE.EXTANT)
    COLL_ROUGH = colless(INTREE.EXTANT.TS)
    message(COLL_ROUGH)
    COLL_ROUGH_SIZESTD = COLL_ROUGH / Ntip(INTREE.EXTANT)
    message(COLL_ROUGH_SIZESTD)
    COLL_HEARD = COLL_ROUGH/((INTREE.EXTANT$Nnode-1)*(INTREE.EXTANT$Nnode-2)/2)
    message(COLL_HEARD)
    COLL_BLUM = colless(INTREE.EXTANT.TS, norm='yule')
    message(COLL_BLUM)
    SACK_BLUM = sackin(INTREE.EXTANT.TS, norm='yule')
    message(SACK_BLUM)
    AB_CF = maxlik.betasplit(INTREE.EXTANT.TS, up=2, confidence.interval='profile', conf.level=0.95, size.bootstrap=10000)
    message(AB_CF$max_lik)
    message(AB_CF$conf_interval[1])
    message(AB_CF$conf_interval[2])
    PHG = gammaStat(INTREE.EXTANT)
    message(PHG)
    
    BTIMES = branching.times(INTREE.EXTANT)
    RTDIST = max(BTIMES)
    message(RTDIST)
    
    RESULTS.TEMP[3] = COLL_ROUGH
    RESULTS.TEMP[4] = COLL_ROUGH_SIZESTD
    RESULTS.TEMP[5] = COLL_HEARD
    RESULTS.TEMP[6] = COLL_BLUM
    RESULTS.TEMP[7] = SACK_BLUM
    RESULTS.TEMP[8] = AB_CF$max_lik
    RESULTS.TEMP[9] = AB_CF$conf_interval[1]
    RESULTS.TEMP[10] = AB_CF$conf_interval[2]
    RESULTS.TEMP[11] = PHG
    RESULTS.TEMP[12] = RTDIST
    
  }
  return(RESULTS.TEMP)
}
#####

#####begin main loop
for (seed in 1:SVSIZE[1])
{
  SEEDSTR = as.character(SEED_VECTOR[seed])
  message(SEEDSTR)
  FILEPREFIX1 = paste('s',SEEDSTR,'_t',sep="")  #put together first part of input file name
  
  #make names of tar archives containing the Nexus files and put them in a char matrix to step through
  TARNAME1=paste('s',SEEDSTR, '_CTRL.tar',sep='')
  TARNAME2=paste('s',SEEDSTR, '_RANDEXT',EXT_STRENGTH_STR,'.tar',sep='')
  TARNAME3=paste('s',SEEDSTR, '_LOTRTEXT',EXT_STRENGTH_STR,'.tar',sep='')
  TARNAME4=paste('s',SEEDSTR, '_HITRTEXT',EXT_STRENGTH_STR,'.tar',sep='')
  INTARNAMES = rbind(TARNAME1,TARNAME2,TARNAME3,TARNAME4)
  
  #make names of tar archives containing the Nexus files and put them in a char matrix to step through
  OUTFILENAME1=paste('MESA_X5Balance_test_s',SEEDSTR,'_noext',EXT_STRENGTH_STR,'_NEW.out', sep="")
  OUTFILENAME2=paste('MESA_X5Balance_test_s',SEEDSTR,'_randext',EXT_STRENGTH_STR,'_NEW.out', sep="")
  OUTFILENAME3=paste('MESA_X5Balance_test_s',SEEDSTR,'_lotrtext',EXT_STRENGTH_STR,'_NEW.out', sep="")
  OUTFILENAME4=paste('MESA_X5Balance_test_s',SEEDSTR,'_hitrtext',EXT_STRENGTH_STR,'_NEW.out', sep="")
  OUTFILENAMES = rbind(OUTFILENAME1,OUTFILENAME2,OUTFILENAME3,OUTFILENAME4)
  
  for (ser in 1:SERIESSIZE)  #for each of the 4 treatment types
  {
    
    #first un-archive the tar file
    message("Un-archiving file...")
    CMD = paste("tar -xvf ", INTARNAMES[ser], sep = " ")
    message(CMD)
    system(CMD, ignore.stdout=TRUE, ignore.stderr=TRUE)
    
    #now unzip the Nexus files
    message("Unzipping Nexus files...")
    CMD = "gunzip *.nex.gz"
    message(CMD)
    system(CMD, ignore.stdout=TRUE, ignore.stderr=TRUE)
    
    #create a matrix that will hold all of the output from each time point
    RESULTS.MATRIX = matrix(data=0,nrow=1,ncol=SVLEN)
    
    for (time in 2:TVSIZE[1])
    {
      #time=4
      FILETIME = as.character(TIME_VECTOR[time])
      INTREEFILENAME = paste(FILEPREFIX1,FILETIME,NEXFILEEXT,sep="")
      message(INTREEFILENAME)
      
      INTREE = GetTreeFromNexusFile(INTREEFILENAME)
      INTREE.EXTANT = DropExtinctLeaves(INTREE)
      ThisTreeBalanceStats = CalcBalanceStats(TIME_VECTOR[time],INTREE.EXTANT)
      RESULTS.MATRIX = rbind(RESULTS.MATRIX,ThisTreeBalanceStats)
    }
    
    #when done with all time points, write the output file
    message(OUTFILENAMES[ser])
    write.table(RESULTS.MATRIX,file=OUTFILENAMES[ser],sep=" ")
    
    #Delete unzipped Nexus files when finished before going to next loop iteration
    message("Deleting Nexus files... ")
    system2('rm',args="*.nex", stdout=NULL, stderr=NULL)  #used "rm" hack for DOS because "del" keeps throwing error code 127
  }
}
#####end main loop