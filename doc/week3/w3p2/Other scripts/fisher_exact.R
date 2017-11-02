DATA = as.matrix(c(6,9,0,0,1,0,5,1,0,0,8,18,32,0,0,0,0,7,1,0,4,0,0,0,1,4,2,0,0,0,0,0,0,0,1))
dim(DATA) = c(5,7)
fisher.test(DATA,workspace = 20000000)

DATA =as.matrix(c(9,5,23,2,4,5,1,4,1,25,6,1,0,0,2,0,8,0,0,1,0,1,0,2,0,0,0,0,0,0,0,0,0,0,0))
dim(DATA) = c(7,5)
fisher.test(DATA,workspace = 20000000)

