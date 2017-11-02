setwd("C:/Users/Gabriel Yedid/Desktop/YGD_MESAPROJECT/NEWOUTPUTFILES/")

SEED_LIST_ALL = read.table("MESA_seed_list.txt",header=F)

#only use to get rid of any missing seed values
#SEED_LIST_ALL = SEED_LIST_ALL[-5,1]
#SEED_LIST_ALL = SEED_LIST_ALL[-8]
#SEED_LIST_ALL = SEED_LIST_ALL[-14]
#SEED_LIST_ALL = as.data.frame(SEED_LIST_ALL)

hndl1 = file("LOTRT_0.5_rootrep_times.out", "w")
hndl2 = file("LOTRT_0.5_rootrep_rootages.out", "w")

for (x in 1:dim(SEED_LIST_ALL)[1]) {
  SEEDVAL = as.character(SEED_LIST_ALL[x,1])
  message(SEEDVAL)
  INFILENAME = paste("MESA_X5Balance_test_s",SEEDVAL,"_lotrtext0.5_NEW.out",sep="")
  message(INFILENAME)
  
  DATA = read.table(INFILENAME, header=T)
  RIDROW = which(DATA[,1]==301)
  DATA = DATA[-RIDROW,]
  
  DATADIM = dim(DATA)
  
  ROOT_AGE_VEC = numeric(DATADIM[1])
  TIME_DIFF_VEC = numeric(DATADIM[1])
  
  ROOT_AGE_VEC[1:length(ROOT_AGE_VEC)] = DATA[1:length(ROOT_AGE_VEC),1] - DATA[1:length(ROOT_AGE_VEC),12] 
  
  for (i in (length(TIME_DIFF_VEC)):2)
  {
    TIME_DIFF_VEC[i] = DATA[i,12] - DATA[i-1,12]
  }
  
  # CB_L2.5 = -1.3221
  # CB_U97.5 = 1.834
  # CB_L25 = -0.5747
  # CB_U75 = 0.4669
  # 
  # plot(DATA[,1],DATA[,6],ylim = c(-3.5,(max(DATA[,6])+1)), type="l")
  # 
  # abline(h=0, untf=FALSE)
  # abline(h=CB_L2.5, untf=FALSE, lty=2)
  # abline(h=CB_U97.5, untf=FALSE, lty=2)
  # abline(h=CB_L25, untf=FALSE, lty=3)
  # abline(h=CB_U75, untf=FALSE, lty=3)
  # 
  BALBREAKS = DATA[which(TIME_DIFF_VEC < 4.9),1]
  BALBREAK_RT_AGS = ROOT_AGE_VEC[which(TIME_DIFF_VEC < 4.9)]
  # 
  # for (i in 1:length(BALBREAKS)){
  #    abline(v=BALBREAKS[i],untf=FALSE,lty=3)
  # }
  
  cat(SEEDVAL,BALBREAKS,file=hndl1,sep="\t")
  cat("\n",file=hndl1)
  
  cat(SEEDVAL,BALBREAK_RT_AGS,file=hndl2,sep="\t")
  cat("\n",file=hndl2)

}

close(hndl1)
close(hndl2)