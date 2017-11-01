library(ape)
library(apTreeshape)
library(phytools)

setwd("C:/R_MESADATA")

#TIME_VECTOR = c(0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200)
#TIME_VECTOR = c(0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150,155,160,165,170,175,180,185,190,195,200)
#TIME_VECTOR = c(201,205,210,215,220,225,230,235,240,245,250,255,260,265,270,275,280,285,290,295,300,305,310,315,320,325,330,335,340,345,350,355,360,365,370,375,380,385,390,395,400)
#TIME_VECTOR = c(201,205,210,215,220,225,230,235,240,245,250,255,260,265,270,275,280,285,290,295,300,305,310,315,320,325,330,335,340,345,350,355,360,365,370,375,380,385,390,395,400,405,410,415,420,425,430,435,440,445,450,455,460,465,470,475,480,485,490,495,500)
#TIME_VECTOR = c(0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,201,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,360,370,380,390,400)
#TIME_VECTOR = c(0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150,155,160,165,170,175,180,185,190,195,200,205,210,215,220,225,230,235,240,245,250,255,260,265,270,275,280,285,290,295,300,305,310,315,320,325,330,335,340,345,350,355,360,365,370,375,380,385,390,395,400,405,410,415,420,425,430,435,440,445,450,455,460,465,470,475,480,485,490,495,500,505,510,515,520,525,530,535,540,545,550,555,560,565,570,575,580,585,590,595,600)

#time vector, 600 units
TIME_VECTOR = c(0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150,155,160,165,170,175,180,185,190,195,200,205,210,215,220,225,230,235,240,245,250,255,260,265,270,275,280,285,290,295,300,301,305,310,315,320,325,330,335,340,345,350,355,360,365,370,375,380,385,390,395,400,405,410,415,420,425,430,435,440,445,450,455,460,465,470,475,480,485,490,495,500,505,510,515,520,525,530,535,540,545,550,555,560,565,570,575,580,585,590,595,600)

#time vector, 1000 units
#TIME_VECTOR = c(0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150,155,160,165,170,175,180,185,190,195,200,205,210,215,220,225,230,235,240,245,250,255,260,265,270,275,280,285,290,295,300,305,310,315,320,325,330,335,340,345,350,355,360,365,370,371,375,380,385,390,395,400,405,410,415,420,425,430,435,440,445,450,455,460,465,470,475,480,485,490,495,500,505,510,515,520,525,530,535,540,545,550,555,560,565,570,575,580,585,590,595,600,605,610,615,620,625,630,635,640,645,650,655,660,665,670,675,680,685,690,695,700,705,710,715,720,725,730,735,740,745,750,755,760,765,770,775,780,785,790,795,800,805,810,815,820,825,830,835,840,845,850,855,860,865,870,875,880,885,890,900,905,910,915,920,925,930,935,940,945,950,955,960,965,970,975,980,985,990,995,1000)


TVSIZE = length(TIME_VECTOR)

#tree size extant
#colless_rough
#colless_heard
#colless_blum
#sackin_blum
#aldous B
#B lower CI
#B upper CI
#PHG

STATSVEC = c("TIME", "extant size", "COLL_ROUGH",
             "COLL_ROUGH_SIZESTD", "COLL_HEARD", "COLL_BLUM", "SACK_BLUM",
             "ALDOUS_B", "B_LCI", "B_UCI", "PHG")

NUMSTATS = 10
RESULTS.ARRAY = (1:(TVSIZE[1]*(NUMSTATS+1)))
dim(RESULTS.ARRAY)=c(TVSIZE[1],(NUMSTATS+1))
RESULTS.ARRAY[1:TVSIZE,1:(NUMSTATS+1)] = 0
colnames(RESULTS.ARRAY) = STATSVEC

SEEDSTR = '4982'
FILEPREFIX1 = paste('s',SEEDSTR,'_t',sep="")
#FILEPREFIX1 = 'out4dgPX_t'
#FILEPREFIX1 = 'out4aaNOEXT_t'
#FILEPREFIX1 = 'out4aEXT_t'
FILEEXT = '_1.nex'
for (time in 2:TVSIZE[1])
  #for (time in 3:TVSIZE[1])
{
  #time = 2
  print(time)
  print(TIME_VECTOR[time])
  RESULTS.ARRAY[time,1] = TIME_VECTOR[time]
  #FILETIME = as.character(TIME_VECTOR[time])
  FILETIME = as.character(TIME_VECTOR[time])
  
  INTREEFILENAME = paste(FILEPREFIX1,FILETIME,FILEEXT,sep="")
  #INTREEFILENAME = ('out4aPXNOEXT_t355_1.nex')
  print(INTREEFILENAME)
  #----------extract the Newick structure from the Nexus file-----#
  INFILE = file(INTREEFILENAME,"r")
  INLINE = readLines(INFILE,n=1)  #get 1st line to make sure not stuck
  #INFILE = file(IFNAME)
  print('Finding tree structure...')
  while(TRUE)
  {
    INLINE = readLines(INFILE,n=1)
    if (length(INLINE) == 0) {
      print('reached end of file')
      break
    }
    else {
      #print(INLINE)
      
      TEST = grep("TREE * ", INLINE)
      if (length(TEST) > 0) {
        GRABLINE = INLINE
        print('Got tree structure!')
      }
    }
  }
  close(INFILE)
  #print(GRABLINE)
  
  GRABSPLIT = strsplit(GRABLINE,split="] ")
  TREESTRUCT = GRABSPLIT[[1]][2]
  #outfilename = 'tmptree.txt'
  #OUTFILE = file(outfilename,"w")
  #writeLines(TREESTRUCT,OUTFILE)
  #close(OUTFILE)
  INTREE = read.tree(text=TREESTRUCT)
  #INTREE = read.nexus(INTREEFILENAME)
  #-------------------------------------------------------------#
  
  #INTREE = read.nexus("s4995_t305_1.nex")
  summary(INTREE)
  #    plot(INTREE, show.tip.label=FALSE)
  #    ltt.plot(INTREE)
  
  
  print('Getting extincts...')
  tips = getExtinct(INTREE, tol=0.005)
  Ntip(INTREE) - length(tips)
  print('Dropping extincts...')
  INTREE.EXTANT = drop.tip(INTREE,tips)
  summary(INTREE.EXTANT)
  #     plot(INTREE.EXTANT, show.tip.label=FALSE)
  #     ltt.plot(INTREE.EXTANT)
  
  
  if (Ntip(INTREE.EXTANT) >= 3)
  {
    print('Calculating tree shape stats...')
    INTREE.EXTANT.TS = as.treeshape(INTREE.EXTANT)
    COLL_ROUGH = colless(INTREE.EXTANT.TS)
    print(COLL_ROUGH)
    COLL_ROUGH_SIZESTD = COLL_ROUGH / Ntip(INTREE.EXTANT)
    print(COLL_ROUGH_SIZESTD)
    COLL_HEARD = COLL_ROUGH/((INTREE.EXTANT$Nnode-1)*(INTREE.EXTANT$Nnode-2)/2)
    print(COLL_HEARD)
    COLL_BLUM = colless(INTREE.EXTANT.TS, norm='yule')
    print(COLL_BLUM)
    #colless.test(INTREE.EXTANT.TS, model="yule",alternative="greater", n.mc=10000)
SACK_BLUM = sackin(INTREE.EXTANT.TS, norm='yule')
print(SACK_BLUM)
AB_CF = maxlik.betasplit(INTREE.EXTANT.TS, up=2,
                         confidence.interval='profile', conf.level=0.95, size.bootstrap=10000)
print(AB_CF$max_lik)
print(AB_CF$conf_interval[1])
print(AB_CF$conf_interval[2])
PHG = gammaStat(INTREE.EXTANT)
print(PHG)
RESULTS.ARRAY[time,2] = Ntip(INTREE.EXTANT)
RESULTS.ARRAY[time,3] = COLL_ROUGH
RESULTS.ARRAY[time,4] = COLL_ROUGH_SIZESTD
RESULTS.ARRAY[time,5] = COLL_HEARD
RESULTS.ARRAY[time,6] = COLL_BLUM
RESULTS.ARRAY[time,7] = SACK_BLUM
RESULTS.ARRAY[time,8] = AB_CF$max_lik
RESULTS.ARRAY[time,9] = AB_CF$conf_interval[1]
RESULTS.ARRAY[time,10] = AB_CF$conf_interval[2]
RESULTS.ARRAY[time,11] = PHG
  }
  
}

plot(INTREE, show.tip.label=FALSE)
ltt.plot(INTREE)
plot(INTREE.EXTANT, show.tip.label=FALSE)
ltt.plot(INTREE.EXTANT)

TVEC = RESULTS.ARRAY[,1]
TSVEC = RESULTS.ARRAY[,2]
CRSTDVEC = RESULTS.ARRAY[,4]
CHVEC = RESULTS.ARRAY[,5]
CBVEC = RESULTS.ARRAY[,6]
SBVEC = RESULTS.ARRAY[,7]
ABVEC = RESULTS.ARRAY[,8]
PHGVEC = RESULTS.ARRAY[,11]

plot(TVEC,TSVEC)

plot(TVEC,CRSTDVEC)

plot(TVEC,log10(CHVEC))

UPPER_YLIM = max(CBVEC)+1
plot(TVEC,CBVEC, ylim=c(-2,UPPER_YLIM))
abline(h=0, untf=FALSE)
abline(h=-1.33, untf=FALSE, lty=2)
abline(h=1.9, untf=FALSE, lty=2)


plot(TVEC,SBVEC)
abline(h=0, untf=FALSE)

plot(TVEC,ABVEC, ylim=c(-2,2))
lines(TVEC,RESULTS.ARRAY[,9],lty=1, lwd=1)
lines(TVEC,RESULTS.ARRAY[,10],lty=1, lwd=1)
abline(h=0, untf=FALSE)

plot(TVEC, PHGVEC)
abline(h=0, untf=FALSE)

#write.table(RESULTS.ARRAY, file ='MESA_X5Balance_test_s4983_randext0.9.out', sep=" ")
#write.table(RESULTS.ARRAY, file ='MESA_X5Balance_test_s4983_lotrext0.9.out', sep=" ")
#write.table(RESULTS.ARRAY, file ='MESA_X5Balance_test_s4982_hitrext0.9MODIFIED.out', sep=" ")
#write.table(RESULTS.ARRAY, file ='MESA_X5Balance_test_s5000_noext0.9.out', sep=" ")

# YULEDIST_SIZE800 = sort(sapply(rtreeshape(25000,tip.number=800,model="yule"),FUN=colless,norm="yule"))
# YULEDIST_SIZE700 = sort(sapply(rtreeshape(25000,tip.number=700,model="yule"),FUN=colless,norm="yule"))
# YULEDIST_SIZE600 = sort(sapply(rtreeshape(25000,tip.number=600,model="yule"),FUN=colless,norm="yule"))
# YULEDIST_SIZE500 = sort(sapply(rtreeshape(25000,tip.number=500,model="yule"),FUN=colless,norm="yule"))
# YULEDIST_SIZE400 = sort(sapply(rtreeshape(25000,tip.number=400,model="yule"),FUN=colless,norm="yule"))
# YULEDIST_SIZE768 = sort(sapply(rtreeshape(25000,tip.number=768,model="yule"),FUN=colless,norm="yule"))
# 
# 
# mean(YULEDIST_SIZE768)
# YULEDIST_SIZE768[25000*0.025]
# YULEDIST_SIZE768[25000*0.975]
# YULEDIST_SIZE768[25000*0.25]
# YULEDIST_SIZE768[25000*0.75]
# median(YULEDIST_SIZE768)
# min(YULEDIST_SIZE768)
# max(YULEDIST_SIZE768)
