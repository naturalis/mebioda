setwd("C:/Users/Gabriel Yedid/Desktop/YGD_MESAPROJECT/NEWOUTPUTFILES/")
#setwd("C:/R_MESADATA")

CB_L2.5 = -1.3221
CB_U97.5 = 1.834
CB_L25 = -0.5747
CB_U75 = 0.4669

#RESULTS.MATRIX.CTRL = read.table("MESA_X5Balance_test_s4983_noex0tBEXT0.05_PROLONG1200_CHANGETRAITLIM.out", header=T)

RESULTS.MATRIX.CTRL = read.table("MESA_X5Balance_test_s4948_noext0.9_NEW.out")

TIME_VEC = RESULTS.MATRIX.CTRL[,1]
CBVEC_CTRL = RESULTS.MATRIX.CTRL[,6]#balance
ROOTAGEVEC =  RESULTS.MATRIX.CTRL[,1] - RESULTS.MATRIX.CTRL[,12]
TIMEANDAGEVEC = cbind(TIME_VEC,ROOTAGEVEC)

ROOTAGEJUMPSTARTS = c(330,335,340,360,500,560)  #for s4952!!!
ROOTAGEJUMPENDS = c(335,340,345,365,505,565)  #for s4952!!!

ROOTAGEJUMPSTARTS = c(440,475,490)  #for s4918!!!
ROOTAGEJUMPENDS = c(445,480,495)  #for s4918!!!

ROOTAGEJUMPSTARTS = c(455,510,520,545,560,580)  #for s4948!!!
ROOTAGEJUMPENDS = c(460,515,525,550,565,585)  #for s4948!!!


#########TO PLOT 2 GRAPHS ONE ON TOP OF THE OTHER
#########LOWER PLOT IS EXTANT SIZE VS. TIME
#########UPPER PLOT HERE IS COLLESS_BLUM VS TIME, X AXIS SUPPRESSED

par(fig=c(0,1,0.5,1))  #define upper plotting region
par(mar=c(0,5,2,2))  #leave off lower margin (1st argument)


UPPER_YLIM = 16
LOWER_YLIM = -2.5
plot(TIME_VEC,CBVEC_CTRL,ylim =c(LOWER_YLIM,UPPER_YLIM),main="",xlab="",xaxt="n",xaxs="i",yaxs="i",ylab="I_c(Blum)",type="n")

abline(h=0, untf=FALSE)
abline(h=CB_L2.5, untf=FALSE, lty=2)
abline(h=CB_U97.5, untf=FALSE, lty=2)
abline(h=CB_L25, untf=FALSE, lty=3)
abline(h=CB_U75, untf=FALSE, lty=3)

points(TIME_VEC,CBVEC_CTRL,type="l",col="black")

#abline(v=c(25,50,75,100,125,150,175,200,225,250,275,325,350,375,400,425,450,475,500,525,550,575,600),untf=FALSE, lty=2, lwd=0.5, col="lightblue")
for (x in 1:length(TIME_VEC))
{
  if ((TIME_VEC[x] %% 25)== 0){
    message(TIME_VEC[x])
    abline(v=TIME_VEC[x],untf=FALSE, lty=2, lwd=0.5, col="lightblue") }
}
abline(v=300,untf=FALSE, lty=4)

for (x in 1:length(ROOTAGEJUMPSTARTS))
{
    abline(v=ROOTAGEJUMPSTARTS[x],untf=FALSE, lty=2, lwd=0.5, col="darkgrey")
}

for (x in 1:length(ROOTAGEJUMPENDS))
{
  abline(v=ROOTAGEJUMPENDS[x],untf=FALSE, lty=2, lwd=0.5, col="black")
}


legend("topleft",c("CTRL"),lty=c(1,1,1,1),col=c("black"), bty="n",ncol=1)
box()


#####now plot lower portion
par(fig=c(0,1,0,0.5),new=TRUE)  #define lower plotting region, adding to current figure
par(mar=c(5,5,0,2))  #leave off upper margin (3rd argument)

# TSVEC_CTRL = RESULTS.MATRIX.CTRL[,2]
# TSVEC_RANDEXT = RESULTS.MATRIX.RANDEXT[,2]
# TSVEC_LOTRTEXT = RESULTS.MATRIX.LOTRTEXT[,2]
# TSVEC_HITRTEXT = RESULTS.MATRIX.HITRTEXT[,2]
# plot(TIME_VEC,TSVEC_CTRL,ylim =c(0,768),main="",xlab="time",ylab="extant size",type="n")
# 
# points(TIME_VEC,TSVEC_RANDEXT,type="l",col="red")
# points(TIME_VEC,TSVEC_LOTRTEXT,type="l",col="blue")
# points(TIME_VEC,TSVEC_HITRTEXT,type="l",col="purple")
# points(TIME_VEC,TSVEC_CTRL,type="l",col="black")
# abline(v=300,untf=FALSE, lty=4)
# abline(v=c(25,50,75,100,125,150,175,200,225,250,275,325,350,375,400,425,450,475,500,525,550,575),untf=FALSE, lty=2, lwd=0.5, col="lightblue")
# 
# legend("bottomright",c("ctrl","randext","lotrtext","hitrtext"),lty=c(1,1,1,1),col=c("black","red","blue","purple"),bty="n",ncol=2)

ROOTAGE_YLIMS = c(0,max(ROOTAGEVEC)+10)

plot(TIME_VEC,ROOTAGEVEC,main="",ylim=ROOTAGE_YLIMS,xlab="time",xaxs="i",yaxs="i",ylab="Root Age",type="n")

points(TIME_VEC,ROOTAGEVEC,type="l",col="black")

#abline(v=c(25,50,75,100,125,150,175,200,225,250,275,325,350,375,400,425,450,475,500,525,550,575,600),untf=FALSE, lty=2, lwd=0.5, col="lightblue")
for (x in 1:length(TIME_VEC))
{
  if ((TIME_VEC[x] %% 25)== 0){
    message(TIME_VEC[x])  
    abline(v=TIME_VEC[x],untf=FALSE, lty=2, lwd=0.5, col="lightblue") }
}
abline(v=300,untf=FALSE, lty=4)

for (x in 1:length(ROOTAGEJUMPSTARTS))
{
  abline(v=ROOTAGEJUMPSTARTS[x],untf=FALSE, lty=2, lwd=0.5, col="darkgrey")
}

for (x in 1:length(ROOTAGEJUMPENDS))
{
  abline(v=ROOTAGEJUMPENDS[x],untf=FALSE, lty=2, lwd=0.5, col="black")
}

abline(h=0, untf=FALSE)
box()

