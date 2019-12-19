

# This tutorial looks at elapid snakes (venomous snakes)
# Does adaptation to the marine environment lead to faster or slower speciation rates?
# i.e. to sea snakes (and sea kraits) speciate more rapidly or more slowly than their terrestrial relatives
# The tutorial uses a phylogenetic tree of elapids, and a table of their habitats
# Binary State Speciation and Extinction (BiSSE) is then used to assess whether marine snakes speciate faster

### load packages
library(diversitree)
library(phangorn)

#read in the tree and habitat data for each taxon
#time-scaled molecular phylogeny of venomous snakes
tree <- read.nexus("elapids.tree") 
#a file which notes taxa that are terrestrial and those that are aquatic
df <- read.csv("habitats.csv") 

#do some stuff to get habitat data in the correct format
v <- df[,2]
names(v) <- df[,1]

# set up the BiSSE model
# This is the full model with 6 parameters
# Each state (terrestrial and marine) has a separate speciation and extinction rate
# There are also parameters for the rate of transition from terrestrial to marine and vice versa

bisse6 <- make.bisse(tree, v)

#get some starting parameters
start <- starting.point.bisse(tree) 

#run the model
ml6 <- find.mle(bisse6, start)

# constrain the BiSSE model
# first a 5 parameter model with equal extinction rates for the two characters
bisse5 <- constrain(bisse6, mu1 ~ mu0) 

# next a 4 parameter model with equal speciation and extinction rates
bisse4 <- constrain(bisse6, lambda1 ~ lambda0, mu1 ~ mu0)


#run the constrained models 
ml5 <- find.mle(bisse5, start)
ml4 <- find.mle(bisse4, start)

#check log likelihood values
ml6$lnLik
ml5$lnLik
ml4$lnLik

#something strange has happened
#How can the likelhood for the 5 parameter model be lower than the constrained?
#appears that the analysis is stuck on a local peak
#try using values from the unconstrained analysis as the starting values for the constrained analysis
start <- ml6$par
ml6.1 <- find.mle(bisse6, start)
ml5.1 <- find.mle(bisse5, start) 
ml4.1 <- find.mle(bisse4, start)

#Has the likelihood value changed?
ml6.1$lnLik
ml5.1$lnLik
ml4.1$lnLik

# Now calculate AIC for each model
2*6 - 2*ml6.1$lnLik
2*5 - 2*ml5.1$lnLik
2*4 - 2*ml4$lnLik

#which model is the best fit?

# There is strong evidence for a difference
# But which way is it?
# Check the parameter values
ml5.1$par


#optional, run ancestral state reconstruction under BiSSE
asr <- asr.marginal(bisse5, coef(ml5.1))
plot(tree, cex=0.4, tip.col=rgb(1-v,0, v), type="fan")
nodelabels(frame="circle", col=rgb(1-asr[2,], 0, asr[2,]), pch=20, cex=1)

# So marine snakes speciate faster according to BiSSE
# Do you believe the result?

############################

#recall that we need to construct appropriate null models
library(hisse)

#
trans.rates.bisse <- TransMatMaker(hidden.states=FALSE)
trans.rates.bisse

#To simplify we will look at equal rates

trans.rates.bisse.eq <- ParEqual(trans.rates.bisse, c(1,2))

#set diversification parameters in order 0A, 1A, 0B, 1B
#This is the BiSSE model

turnover.anc.bisse <- c(1,2,0,0)

# To simplifiy, we run the models without extinction
eps.anc <- c(0,0,0,0)

#run the bisse model
pp.bisse <- hisse(tree, df, f=c(0.42,0.68), hidden.states=FALSE, turnover.anc=turnover.anc.bisse, eps.anc=eps.anc, trans.rate=trans.rates.bisse.eq, output.type="net.div")

#Make an ancestral state reconstruction
pp.bisse.recon <- MarginRecon(phy=tree, data=df, f = pp.bisse$f, pars = pp.bisse$solution, aic = pp.bisse$AIC, n.cores=1)

#plot it
plot(pp.bisse.recon, fsize=0.5)

#make a transition rate matrix with hidden states
trans.rates <-TransMatMaker(hidden.states=TRUE)

#this has the maximum number of possible parameters
#simplify by removing the possibility of simultaneous jumps between hidden states (e.g. 0A -> 0B)

trans.rates.nodual =ParDrop(trans.rates,c(3,5,8,10))

#set all remaining rates to be equal
trans.rates.nodual.allequal <-ParEqual(trans.rates.nodual,c(1,2,1,3,1,4,1,5,1,6,1,7,1,8))
trans.rates.nodual.allequal

#set diversification parameters in order 0A, 1A, 0B, 1B
#This is the null model i.e. state independent diversification
turnover.anc.null <- c(1,1,2,2)

# Run the null model
pp.null <- hisse(tree, df, f=c(0.42,0.68), hidden.states=TRUE, turnover.anc=turnover.anc.null, eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal, output.type="net.div")

#ancestral reconstruction
pp.null.recon <- MarginRecon(phy=tree, data=df, f = pp.null$f, pars = pp.null$solution, aic = pp.null$AIC, n.cores=1)

plot(pp.null.recon, fsize=0.5)

#looks wrong, try new starting values
starting.vals <- c(0.01, 0, 0.002)
pp.null <- hisse(tree, df, f=c(0.42,0.68), hidden.states=TRUE, turnover.anc=turnover.anc.null, eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal, output.type="net.div", starting.vals=starting.vals)
pp.null.recon <- MarginRecon(phy=tree, data=df, f = pp.null$f, pars = pp.null$solution, aic = pp.null$AIC, n.cores=1)

plot(pp.null.recon, fsize=0.5)


hisse.results.list <- list()
hisse.results.list[[1]] <- pp.bisse.recon
hisse.results.list[[2]] <- pp.null.recon
plot.hisse.states(hisse.results.list, rate.param="net.div", fsize=0.5)


############################

# Now let us try the non-parametric method FiSSE

# First load the functions
source("traitDependent_functions.R")

# Run the test
res <- FISSE.binary(tree, v)

# Calculate p value
pval   <- min(res$pval, 1-res$pval)*2
pval

# Does this method give a different result from BiSSE?
# Compare the speciation rates between BiSSE and FiSSE. Are they different?

############################

