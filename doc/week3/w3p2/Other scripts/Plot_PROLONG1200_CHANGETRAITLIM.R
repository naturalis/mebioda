setwd("C:/Users/Gabriel Yedid/Desktop/YGD_MESAPROJECT/NEWOUTPUTFILES/")

#Zeileis et al.'s function for drawing palettes
pal <- function(col, border = "light gray", ...)
{
  n <- length(col)
  plot(0, 0, type="n", xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "", ...)
  rect(0:(n-1)/n, 0, 1:n/n, 1, col = col, border = border)
}



DATA2 = read.table("noext0BEXT0.05_PROLONG1200_CHANGETRAITLIM_BalStats_COLBLUM_Summary.out",header=T,sep=" ")
DATA1 = read.table("noext0BEXT0.05_PROLONG1200_CHANGETRAITLIM_TraitVar_SummaryMOD.out.txt",header=T,sep="\t")

DATA1 = DATA1[,-1]
colnames(DATA1) = colnames(DATA2)

PAL = rainbow(30)

REPNAMES = colnames(DATA1)
REPNAMES = REPNAMES[-1]

MAXBALVALS1 = vector(mode="numeric",length=length(REPNAMES))
MAXVARVALS1 = vector(mode="numeric",length=length(REPNAMES))
MAXBALTIMES1 = vector(mode="numeric",length=length(REPNAMES))
MAXVARTIMES1 = vector(mode="numeric",length=length(REPNAMES))

MAXBALVALS2 = vector(mode="numeric",length=length(REPNAMES))
MAXVARVALS2 = vector(mode="numeric",length=length(REPNAMES))
MAXBALTIMES2 = vector(mode="numeric",length=length(REPNAMES))
MAXVARTIMES2 = vector(mode="numeric",length=length(REPNAMES))

TSWITCHLIM = 161

for (n in 1:length(REPNAMES)) {
  
  #n=1
  MXBV1 = max(DATA2[1:TSWITCHLIM,n+1])
  MXBV1ROW = which(DATA2[1:TSWITCHLIM,n+1]==MXBV1)
  MXBV1TIME = DATA1[MXBV1ROW,1]
  
  MXTV1 = max(DATA1[1:TSWITCHLIM,n+1])
  MXTV1ROW = which(DATA1[1:TSWITCHLIM,n+1]==MXTV1)
  MXTV1TIME = DATA1[MXTV1ROW,1]
  
  MAXBALVALS1[n] = MXBV1
  MAXBALTIMES1[n] = MXBV1TIME
  MAXVARVALS1[n] = MXTV1
  MAXVARTIMES1[n] = MXTV1TIME
  
  MXBV2 = max(DATA2[TSWITCHLIM:241,n+1])
  MXBV2ROW = which(DATA2[TSWITCHLIM:241,n+1]==MXBV2) + (TSWITCHLIM-1)
  MXBV2TIME = DATA1[MXBV2ROW,1]
  
  MXTV2 = max(DATA1[TSWITCHLIM:241,n+1])
  MXTV2ROW = which(DATA1[TSWITCHLIM:241,n+1]==MXTV2) + (TSWITCHLIM-1)
  MXTV2TIME = DATA1[MXTV2ROW,1]
  
  MAXBALVALS2[n] = MXBV2
  MAXBALTIMES2[n] = MXBV2TIME
  MAXVARVALS2[n] = MXTV2
  MAXVARTIMES2[n] = MXTV2TIME

} 

MAXFRAME = data.frame(REPNAMES,MAXBALVALS1,MAXBALTIMES1,MAXVARVALS1,MAXVARTIMES1,MAXBALVALS2,MAXBALTIMES2,MAXVARVALS2,MAXVARTIMES2)


#####now plot upper portion
par(fig=c(0,1,0.5,1))  #define upper plotting region
par(mar=c(0,5,2,2))  #leave off lower margin (1st argument)

VAR_VEC = DATA1[,2]

TIME_VEC = DATA1[,1]
plot(TIME_VEC,VAR_VEC,ylim = c(0,5),main="",xaxs="i",xaxt="n",xlab="time",ylab="trait variance",type="n")
points(TIME_VEC,VAR_VEC,type="l",col=PAL[1])

for (n in 1:length(PAL))
{
  VAR_VEC = DATA1[,n+1]
  
  #plot(TIME_VEC,VAR_VEC,ylim = c(0,6),main="",xaxs="i",yaxs="i",xaxt="n",xlab="time",ylab="trait variance",type="n")
  points(TIME_VEC,VAR_VEC,type="l",col=PAL[n])
  abline(v=MAXFRAME[n,5],untf=FALSE,lty=2,col=PAL[n])
  abline(v=MAXFRAME[n,9],untf=FALSE,lty=2,col=PAL[n])
  
  
}


par(fig=c(0,1,0,0.5),new=TRUE)  #define lower plotting region, adding to current figure
par(mar=c(5,5,0,2))  #leave off upper margin (3rd argument)

UPPER_YLIM = 15
LOWER_YLIM = -2.5

CB_L2.5 = -1.3221
CB_U97.5 = 1.834
CB_L25 = -0.5747
CB_U75 = 0.4669


BAL_VEC = DATA2[,2]

plot(TIME_VEC,BAL_VEC,ylim = c(LOWER_YLIM,UPPER_YLIM),main="",xaxs="i",yaxs="r",xlab="time",ylab="I_c(Blum)",type="n")
points(TIME_VEC,BAL_VEC,type="l",col=PAL[1])

for (n in 1:length(PAL))
{
  BAL_VEC = DATA2[,n+1]
  
  #plot(TIME_VEC,VAR_VEC,ylim = c(0,6),main="",xaxs="i",yaxs="i",xaxt="n",xlab="time",ylab="trait variance",type="n")
  points(TIME_VEC,BAL_VEC,type="l",col=PAL[n])
  abline(v=MAXFRAME[n,3],untf=FALSE,lty=2,col=PAL[n])
  abline(v=MAXFRAME[n,7],untf=FALSE,lty=2,col=PAL[n])
  
  
}



abline(h=0, untf=FALSE)
abline(h=CB_L2.5, untf=FALSE, lty=2)
abline(h=CB_U97.5, untf=FALSE, lty=2)
abline(h=CB_L25, untf=FALSE, lty=3)
abline(h=CB_U75, untf=FALSE, lty=3)










