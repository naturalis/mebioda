setwd("C:/CUSTOM_R_FILES")
#setwd("C:/R_MESADATA")

CB_L2.5 = -1.3221
CB_U97.5 = 1.834
CB_L25 = -0.5747
CB_U75 = 0.4669

#RESULTS.MATRIX.CTRL = read.table("MESA_X5Balance_test_s4983_noex0tBEXT0.05_PROLONG1200_CHANGETRAITLIM.out", header=T)

RESULTS.MATRIX.CTRL = read.table("MESA_X5Balance_test_s4997_noext0.9.out")
RESULTS.MATRIX.RANDEXT = read.table("MESA_X5Balance_test_s4997_randext0.9.out")
RESULTS.MATRIX.LOTRTEXT = read.table("MESA_X5Balance_test_s4997_lotrtext0.9.out")
RESULTS.MATRIX.HITRTEXT = read.table("MESA_X5Balance_test_s4997_hitrtext0.9.out")


###############COLLESS_BLUM
UPPER_YLIM = max(max(RESULTS.MATRIX.CTRL[,6]),max(RESULTS.MATRIX.RANDEXT[,6]),max(RESULTS.MATRIX.LOTRTEXT[,6]),max(RESULTS.MATRIX.HITRTEXT[,6]))
LOWER_YLIM = -2.5 #min(min(RESULTS.MATRIX.CTRL[,6]),min(RESULTS.MATRIX.RANDEXT[,6]),min(RESULTS.MATRIX.LOTRTEXT[,6]),min(RESULTS.MATRIX.HITRTEXT[,6]))
TIME_VEC = RESULTS.MATRIX.CTRL[,1]
CBVEC_CTRL = RESULTS.MATRIX.CTRL[,6]
CBVEC_RANDEXT = RESULTS.MATRIX.RANDEXT[,6]
CBVEC_LOTRTEXT = RESULTS.MATRIX.LOTRTEXT[,6]
CBVEC_HITRTEXT = RESULTS.MATRIX.HITRTEXT[,6]
YLAB_EXPR = expression(paste("I"[C],"_Blum"))
#plot(TIME_VEC,CBVEC_CTRL,ylim =c(LOWER_YLIM,UPPER_YLIM),main="COLLESS_BLUM",xlab="time",ylab=YLAB_EXPR,type="n",family="sans",cex.axis=1.25,cex.lab=1.25,yaxs="i",xaxs="i")
plot(TIME_VEC,CBVEC_CTRL,ylim =c(LOWER_YLIM,UPPER_YLIM),xlab="time",ylab=YLAB_EXPR,type="n",family="sans",cex.axis=1.25,cex.lab=1.25,yaxs="i",xaxs="i")

abline(h=0, untf=FALSE)
abline(h=CB_L2.5, untf=FALSE, lty=2)
abline(h=CB_U97.5, untf=FALSE, lty=2)
abline(h=CB_L25, untf=FALSE, lty=3)
abline(h=CB_U75, untf=FALSE, lty=3)
abline(v=300,untf=FALSE, lty=4)

points(TIME_VEC,CBVEC_RANDEXT,type="l",col="red")
points(TIME_VEC,CBVEC_LOTRTEXT,type="l",col="blue")
points(TIME_VEC,CBVEC_HITRTEXT,type="l",col="purple")
points(TIME_VEC,CBVEC_CTRL,type="l",col="black")


#randext=RANDOM     #lotrtext=SOD     #hitrtext=SOR

#for s4918
legend("topleft",c("CTRL","RANDOM","SOD","SOR"),lty=c(1,1,1,1),col=c("black","red","blue","purple"), bty="n",ncol=2)
abline(v=320,col="red",lty=3) #CSR rand
abline(v=315,col="purple",lty=3) #CSR SOR
abline(v=350,col="blue",lty=3) #CSR SOD

#for s4952
legend("topleft",c("CTRL","RANDOM","SOD","SOR"),lty=c(1,1,1,1),col=c("black","red","blue","purple"), bty="n",ncol=1)
abline(v=310,col="red",lty=3) #CSR rand
abline(v=310,col="purple",lty=3) #CSR SOR
abline(v=335,col="blue",lty=3) #CSR SOD

#for s4997
legend("topleft",c("CTRL","RANDOM","SOD","SOR"),lty=c(1,1,1,1),col=c("black","red","blue","purple"), bty="n",ncol=2)
abline(v=345,col="red",lty=3) #CSR rand
abline(v=335,col="purple",lty=3) #CSR SOR
abline(v=360,col="blue",lty=3) #CSR SOD




#**************

###############EXTANT_SIZE
TSVEC_CTRL = RESULTS.MATRIX.CTRL[,2]
TSVEC_RANDEXT = RESULTS.MATRIX.RANDEXT[,2]
TSVEC_LOTRTEXT = RESULTS.MATRIX.LOTRTEXT[,2]
TSVEC_HITRTEXT = RESULTS.MATRIX.HITRTEXT[,2]
plot(TIME_VEC,TSVEC_CTRL,ylim =c(0,768),main="EXTANT_SIZE",xlab="time",ylab="extant size",type="n")

points(TIME_VEC,TSVEC_RANDEXT,type="l",col="red")
points(TIME_VEC,TSVEC_LOTRTEXT,type="l",col="blue")
points(TIME_VEC,TSVEC_HITRTEXT,type="l",col="purple")
points(TIME_VEC,TSVEC_CTRL,type="l",col="black")
abline(v=300,untf=FALSE, lty=4)

legend("bottomright",c("ctrl","randext","lotrtext","hitrtext"),lty=c(1,1,1,1),col=c("black","red","blue","purple"),bty="n",ncol=2)

#*************

###############PHG
PHGVEC_CTRL = RESULTS.MATRIX.CTRL[,11]
PHGVEC_RANDEXT = RESULTS.MATRIX.RANDEXT[,11]
PHGVEC_LOTRTEXT = RESULTS.MATRIX.LOTRTEXT[,11]
PHGVEC_HITRTEXT = RESULTS.MATRIX.HITRTEXT[,11]

GLOBAL_PHG_MAX = max(c(max(PHGVEC_CTRL),max(PHGVEC_RANDEXT),max(PHGVEC_LOTRTEXT),max(PHGVEC_HITRTEXT)))
GLOBAL_PHG_MIN = min(c(min(PHGVEC_CTRL),min(PHGVEC_RANDEXT),min(PHGVEC_LOTRTEXT),min(PHGVEC_HITRTEXT)))

GPMAX_CEIL = ceiling(GLOBAL_PHG_MAX)
GPMIN_FLOOR = floor(GLOBAL_PHG_MIN)

plot(TIME_VEC,PHGVEC_CTRL,ylim =c(GPMIN_FLOOR,GPMAX_CEIL),main="PHG",xlab="time",ylab="PHG",type="n")
points(TIME_VEC,PHGVEC_RANDEXT,type="l",col="red")
points(TIME_VEC,PHGVEC_LOTRTEXT,type="l",col="blue")
points(TIME_VEC,PHGVEC_HITRTEXT,type="l",col="purple")
points(TIME_VEC,PHGVEC_CTRL,type="l",col="black")
abline(v=300,untf=FALSE, lty=4)
abline(h=0,untf=FALSE, lty=1)

legend("topleft",c("ctrl","randext","lotrtext","hitrtext"),lty=c(1,1,1,1),col=c("black","red","blue","purple"),bty="n",ncol=2)

#*************ALDOUS_B
ABVEC_CTRL = RESULTS.MATRIX.CTRL[,8]
ABVEC_RANDEXT = RESULTS.MATRIX.RANDEXT[,8]
ABVEC_LOTRTEXT = RESULTS.MATRIX.LOTRTEXT[,8]
ABVEC_HITRTEXT = RESULTS.MATRIX.HITRTEXT[,8]

ABVECL_CTRL = RESULTS.MATRIX.CTRL[,9]
ABVECL_RANDEXT = RESULTS.MATRIX.RANDEXT[,9]
ABVECL_LOTRTEXT = RESULTS.MATRIX.LOTRTEXT[,9]
ABVECL_HITRTEXT = RESULTS.MATRIX.HITRTEXT[,9]

ABVECU_CTRL = RESULTS.MATRIX.CTRL[,10]
ABVECU_RANDEXT = RESULTS.MATRIX.RANDEXT[,10]
ABVECU_LOTRTEXT = RESULTS.MATRIX.LOTRTEXT[,10]
ABVECU_HITRTEXT = RESULTS.MATRIX.HITRTEXT[,10]

plot(TIME_VEC,ABVEC_CTRL,ylim =c(-2,2),main="ALDOUS'_B",xlab="time",ylab="Aldous_b",type="n")
abline(h=0,untf=FALSE,lty=1)
abline(v=300,untf=FALSE, lty=4)


points(TIME_VEC,ABVEC_RANDEXT,type="l",col="red",lty=1,lwd=2)
points(TIME_VEC,ABVECL_RANDEXT,type="l",col="red",lty=2)
points(TIME_VEC,ABVECU_RANDEXT,type="l",col="red",lty=2)

points(TIME_VEC,ABVEC_LOTRTEXT,type="l",col="blue",lty=1,lwd=2)
points(TIME_VEC,ABVECL_LOTRTEXT,type="l",col="blue",lty=2)
points(TIME_VEC,ABVECU_LOTRTEXT,type="l",col="blue",lty=2)

points(TIME_VEC,ABVEC_HITRTEXT,type="l",col="purple",lty=1,lwd=2)
points(TIME_VEC,ABVECL_HITRTEXT,type="l",col="purple",lty=2)
points(TIME_VEC,ABVECU_HITRTEXT,type="l",col="purple",lty=2)

points(TIME_VEC,ABVEC_CTRL,type="l",col="black",lty=1,lwd=2)
points(TIME_VEC,ABVECL_CTRL,type="l",col="black",lty=2)
points(TIME_VEC,ABVECU_CTRL,type="l",col="black",lty=2)

legend("bottomright",c("ctrl","randext","lotrtext","hitrtext"),lty=c(1,1,1,1),col=c("black","red","blue","purple"),bty="n",ncol=2,xjust=0.5)

#########TO PLOT 2 GRAPHS ONE ON TOP OF THE OTHER
#########LOWER PLOT IS EXTANT SIZE VS. TIME
#########UPPER PLOT HERE IS COLLESS_BLUM VS TIME, X AXIS SUPPRESSED

par(fig=c(0,1,0.5,1))  #define upper plotting region
par(mar=c(0,5,2,2))  #leave off lower margin (1st argument)


UPPER_YLIM = 16
LOWER_YLIM = -2.5
TIME_VEC = RESULTS.MATRIX.CTRL[,1]
CBVEC_CTRL = RESULTS.MATRIX.CTRL[,6]
CBVEC_RANDEXT = RESULTS.MATRIX.RANDEXT[,6]
CBVEC_LOTRTEXT = RESULTS.MATRIX.LOTRTEXT[,6]
CBVEC_HITRTEXT = RESULTS.MATRIX.HITRTEXT[,6]
plot(TIME_VEC,CBVEC_CTRL,ylim =c(LOWER_YLIM,UPPER_YLIM),main="PHG",xlab="",xaxt="n",xaxs="i",yaxs="i",ylab="I_c(Blum)",type="n")

abline(h=0, untf=FALSE)
abline(h=CB_L2.5, untf=FALSE, lty=2)
abline(h=CB_U97.5, untf=FALSE, lty=2)
abline(h=CB_L25, untf=FALSE, lty=3)
abline(h=CB_U75, untf=FALSE, lty=3)

points(TIME_VEC,CBVEC_RANDEXT,type="l",col="red")
points(TIME_VEC,CBVEC_LOTRTEXT,type="l",col="blue")
points(TIME_VEC,CBVEC_HITRTEXT,type="l",col="purple")
points(TIME_VEC,CBVEC_CTRL,type="l",col="black")

#abline(v=c(25,50,75,100,125,150,175,200,225,250,275,325,350,375,400,425,450,475,500,525,550,575,600),untf=FALSE, lty=2, lwd=0.5, col="lightblue")
for (x in 1:length(TIME_VEC))
{
  if ((TIME_VEC[x] %% 25)== 0){
    message(TIME_VEC[x])
    abline(v=TIME_VEC[x],untf=FALSE, lty=2, lwd=0.5, col="lightblue") }
}
abline(v=300,untf=FALSE, lty=4)

legend("topleft",c("ctrl","randext","lotrtext","hitrtext"),lty=c(1,1,1,1),col=c("black","red","blue","purple"), bty="n",ncol=2)
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

PHG_YLIMS = c(GLOBAL_PHG_MIN,GLOBAL_PHG_MAX)
PHGVEC_CTRL = RESULTS.MATRIX.CTRL[,11]
PHGVEC_RANDEXT = RESULTS.MATRIX.RANDEXT[,11]
PHGVEC_LOTRTEXT = RESULTS.MATRIX.LOTRTEXT[,11]
PHGVEC_HITRTEXT = RESULTS.MATRIX.HITRTEXT[,11]
plot(TIME_VEC,PHGVEC_CTRL,ylim = PHG_YLIMS,main="",xlab="time",xaxs="i",yaxs="i",ylab="PHG",type="n")

points(TIME_VEC,PHGVEC_RANDEXT,type="l",col="red")
points(TIME_VEC,PHGVEC_LOTRTEXT,type="l",col="blue")
points(TIME_VEC,PHGVEC_HITRTEXT,type="l",col="purple")
points(TIME_VEC,PHGVEC_CTRL,type="l",col="black")

#abline(v=c(25,50,75,100,125,150,175,200,225,250,275,325,350,375,400,425,450,475,500,525,550,575,600),untf=FALSE, lty=2, lwd=0.5, col="lightblue")
for (x in 1:length(TIME_VEC))
{
  if ((TIME_VEC[x] %% 25)== 0){
  message(TIME_VEC[x])  
  abline(v=TIME_VEC[x],untf=FALSE, lty=2, lwd=0.5, col="lightblue") }
}
abline(v=300,untf=FALSE, lty=4)

abline(h=0, untf=FALSE)
box()

#legend("bottomright",c("ctrl","randext","lotrtext","hitrtext"),lty=c(1,1,1,1),col=c("black","red","blue","purple"),bty="n",ncol=2)
