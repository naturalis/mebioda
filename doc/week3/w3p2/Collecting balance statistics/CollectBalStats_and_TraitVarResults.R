setwd("C:/R_MESADATA/CHANGETRAITLIM_BALSTATS/")
NUMREPS = 30
NUMOBS = 241

RESTAB = matrix(data=0.0,nrow = NUMOBS, ncol = NUMREPS+1)
SEED_LIST = c(4962,4963,4964,4965,4966,4967,4968,4969,4970,4971,4972,4973,4974,4976,4977,4978,4979,4983,4984,4988,4989,4991,4992,4993,4994,4995,4997,4998,4999,5000)
SEEDLSCHR = c("TIME",as.character(SEED_LIST))
colnames(RESTAB) = SEEDLSCHR


SEED = SEED_LIST[1]
SEEDSTR = as.character(SEED)
INFILENAME = paste("MESA_X5Balance_test_s",SEEDSTR,"_wnoext0BEXT0.05_PROLONG1200_CHANGETRAITLIM_BalStats.out",sep="")
DATA = read.table(INFILENAME,header=T,sep=" ")

RESTAB[,1] = DATA[,1]
RESTAB[,2] = DATA[,6]

for (n in 2:length(SEED_LIST)) {
  
  SEED = SEED_LIST[n]
  message(SEED)
  SEEDSTR = as.character(SEED)
  INFILENAME = paste("MESA_X5Balance_test_s",SEEDSTR,"_wnoext0BEXT0.05_PROLONG1200_CHANGETRAITLIM_BalStats.out",sep="")
  DATA = read.table(INFILENAME,header=T,sep=" ")
  
  RESTAB[,(n+1)] = DATA[,6]
  
}

write.table(RESTAB, "noext0BEXT0.05_PROLONG1200_CHANGETRAITLIM_BalStats_COLBLUM_Summary.out", sep = " ")

setwd("D:/MESADATA_STORAGE/PROLONG1200_CHANGETRAITLIM/")


RESTAB = matrix(data=0.0,nrow = NUMOBS, ncol = NUMREPS+1)
colnames(RESTAB) = SEEDLSCHR


for (n in 1:length(SEED_LIST)) {

  SEED = SEED_LIST[n]
  SEEDSTR = as.character(SEED)
  ARCFILENAME = paste("s",SEEDSTR,"_wnoext0BEXT0.05_PROLONG1200_CHANGETRAITLIM_TraitVarResults.tar",sep="")
  ZIPARCFILENAME = paste(ARCFILENAME,".gz",sep="")
  XTFILENAME = paste("s",SEEDSTR,"_TraitVar_results_CHANGETRAITLIM.csv",sep="")
  
  CMD = paste("gunzip ", ZIPARCFILENAME, sep = "")
  message(CMD)
  system(CMD)
  CMD = paste("tar -xvf ", ARCFILENAME, " ", XTFILENAME, sep = "")
  message(CMD)
  system(CMD)
  CMD = paste("gzip ", ARCFILENAME, sep = "")
  message(CMD)
  system(CMD)
  
  DATA = read.csv(XTFILENAME,header=T)
  if (n==1) {RESTAB[,1] = DATA[,1]}
  
  RESTAB[,n+1] = DATA[,2]

  
}
  
  
CMD = "rm *.csv"
system(CMD)

write.table(RESTAB,"noext0BEXT0.05_PROLONG1200_CHANGETRAITLIM_TraitVar_Summary.out",sep = " ")




