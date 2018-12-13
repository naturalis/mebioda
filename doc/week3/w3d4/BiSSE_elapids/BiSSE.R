

# This tutorial looks at elapid snakes (venomous snakes)
# Does adaptation to the marine environment lead to faster or slower speciation rates?
# i.e. to sea snakes (and sea kraits) speciate more rapidly or more slowly than their terrestrial relatives
# The tutorial uses a phylogenetic tree of elapids, and a table of their habitats
# Binary State Speciation and Extinction (BiSSE) is then used to assess whether marine snakes speciate faster

### load packages
install.packages("diversitree")
install.packages("phangorn")
library(diversitree)
library(phangorn)

#read in the tree and habitat data for each taxon
#time-scaled molecular phylogeny of venomous snakes
read.nexus("elapids.tree") -> tree
#a file which notes taxa that are terrestrial and those that are aquatic
read.csv("habitats.csv") -> data

#do some stuff to get habitat data in the correct format
data[,2] -> chars
data[,1] -> names(chars)

# set up the BiSSE model
# This is the full model with 6 parameters
# Each state (terrestrial and marine) has a separate speciation and extinction rate
# There are also parameters for the rate of transition from terrestrial to marine and vice versa

make.bisse(tree, chars) -> bisse6

#get some starting parameters
starting.point.bisse(tree) -> start

#run the model
find.mle(bisse6, start) -> ml6

# constrain the BiSSE model
# first a 5 parameter model with equal extinction rates for the two characters
constrain(bisse6, mu1 ~ mu0) -> bisse5

# next a 4 parameter model with equal speciation and extinction rates
constrain(bisse6, lambda1 ~ lambda0, mu1 ~ mu0) -> bisse4


#run the constrained models 
find.mle(bisse5, start) -> ml5
find.mle(bisse4, start) -> ml4

#check log likelihood values
ml6$lnLik
ml5$lnLik
ml4$lnLik

#something strange has happened
#How can the likelhood for the 5 parameter model be lower than the constrained?
#appears that the analysis is stuck on a local peak
#try using values from the unconstrained analysis as the starting values for the constrained analysis
ml6$par -> start
find.mle(bisse6, start) -> ml6.1
find.mle(bisse5, start) -> ml5.1
find.mle(bisse4, start) -> ml4.1

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

#run ancestral state reconstruction under BiSSE
asr.marginal(bisse6, coef(ml6)) -> asr
plot(tree, cex=0.4, tip.col=rgb(1-chars,0, chars), type="fan")
nodelabels(frame="circle", col=rgb(1-asr[2,], 0, asr[2,]), pch=20, cex=1)

# So marine snakes speciate faster according to BiSSE
# Do you believe the result?

############################

# Now let us try a more recent method called FiSSE

# First load the functions
source("traitDependent_functions.R")

# Run the test
res <- FISSE.binary(tree, chars)

# Calculate p value
pval   <- min(res$pval, 1-res$pval)*2
pval

# Does this method give a different result from BiSSE?
# Compare the speciation rates between BiSSE and FiSSE. Are they different?

############################

