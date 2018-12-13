

# This function checks trees to see if they pass ape ultrametricity test.
# If not, it computes the differential root-to-tip distance across all tips.
# It adds the appropriate quantity to each terminal branch length to ensure that 
# tree passes ultrametric test.
# Note: this is only a valid method of making trees ultrametric when the 
# 	non-ultrametricity is due to small numerical discrepancies, e.g., 
#   rounding or other floating point issues during phylogeny construction.
# 

check_and_fix_ultrametric <- function(phy){
	
	if (!is.ultrametric(phy)){
		
		vv <- vcv.phylo(phy)
		dx <- diag(vv)
		mxx <- max(dx) - dx
		for (i in 1:length(mxx)){
			phy$edge.length[phy$edge[,2] == i] <- phy$edge.length[phy$edge[,2] == i] + mxx[i]
		}
		if (!is.ultrametric(phy)){
			stop("Ultrametric fix failed\n")
		}	
	}
	
	return(phy)
}



logit <- function(x, min=0, max=1){
	p <- (x-min)/(max-min)
	log(p/(1-p))
}

invlogit <- function(x, min=0, max=1)
{
	p <- exp(x)/(1+exp(x))
	p <- ifelse( is.na(p) & !is.na(x), 1, p ) # fix problems with +Inf
	p * (max-min) + min
}

# Gets a vector of initial parameters 
#	for a bisse, birth-death, or mk2 model
#	using the likelihood functions created 
#	in diversitree with make.mk2, make.bisse, or make.bd
# These can be plugged directly into the corresponding
#	likelihood function. However, they are not guaranteed
#	to generate finite log-likelihoods. 
# Arguments:
#	fx: the diversitree likelihood function
#	lmin: the minimum value across all parameters
#	lmax: the maximum value across all parameters

getStartingParamsDiversitree <- function(fx, lmin, lmax){
		
		lamset <- runif(3, lmin, lmax)
		names(lamset) <- c('lambda', paste('lambda', 0:1, sep=''))
		muset <- runif(3, 0, 1) * lamset
		names(muset) <- c('mu', paste('mu', 0:1, sep=''))
		qset <- runif(4, lmin, lmax * 0.2)
		names(qset) <- c('q01', 'q10', 'q12', 'q21')
		parvec <- c(lamset, muset, qset)
	 
		if (length(setdiff(argnames(fx), names(parvec))) > 0){
			stop("Invalid argnames from function\n")
		}
		
		parset <- intersect(names(parvec), argnames(fx))
				
		return(parvec[parset])
}


# A general purpose optimization function
# that optimizes parameters of a diversitree likelihood function.
# The likelihood function must correspond to one of the following models:
#	a) BiSSE (or any constrained submodel)
#	b) birth-death
#	c) mk2 (2 state character only model)
fitDiversitree <- function(fx, nopt=1, lmin = 0.0001, lmax=20.0, MAXBAD = 1000, initscale = 0.1){

 
	for (i in 1:nopt){
		
		badcount <- 0
		
		iv <- getStartingParamsDiversitree(fx, lmin=lmin, lmax=lmax*initscale)
		
		resx <- try(optim(iv ,fx, method='L-BFGS-B', control=list(maxit=1000, fnscale=-1), lower=lmin, upper=lmax), silent=T)
		while (class(resx) == 'try-error'){
		iv <- getStartingParamsDiversitree(fx, lmin=lmin, lmax=lmax*initscale)
			resx <- try(optim(iv , fx, method='L-BFGS-B', control=list(maxit=1000, fnscale=-1), lower=lmin, upper=lmax), silent=T)
			
			badcount <- badcount + 1
			if (badcount > MAXBAD){
				stop("Too many fails in fitDiversitree\n")
			}
		}
		
		if (i == 1){
			best <- resx
		}else{
			if (best$value < resx$value){
				best <- resx
			}
		}
		
	}
	
	fres <- list(pars=best$par, loglik=best$value)
	fres$AIC <- -2*fres$loglik + 2*length(argnames(fx))
	fres$counts <- best$counts
	#fres$like_function <- fx
	fres$convergence <- best$convergence
	fres$message <- best$message
	return(fres)
}


fitDiversitree_allmodels <- function(tree, traits, index=1, nopt=5){
	# lam0, lam1, mu0, mu1, q01, q10 
	# lam0, lam1, mu0,      q01, q10
	# lam0, lam1, mu0       q01, q10
	# lam0,       mu0,      q01, q10
	
	
	
	
	lfx6 <- make.bisse(tree, traits)
	lfx5lv <- constrain(lfx6, formulae=list(mu0 ~ mu1))
	lfx5mv <- constrain(lfx6, formulae = list(lambda0 ~ lambda1))
	lfx4 <- constrain(lfx5lv, formulae= list(lambda0 ~ lambda1))
	
	# assume max lambda is 30x pb rate	
	lpb <- 30* (length(tree$tip.label) - 2 ) / sum(tree$edge.length)
 
	fitx6 <- try(fitDiversitree(lfx6, nopt=nopt, lmax=lpb))
	fitx5lv <- try(fitDiversitree(lfx5lv, nopt=nopt, lmax=lpb))
	fitx5mv <- try(fitDiversitree(lfx5mv, nopt = nopt, lmax=lpb))
	fitx4 <- try(fitDiversitree(lfx4, nopt = nopt, lmax = lpb))
	  
	resmat <- matrix(NA, nrow=4, ncol=12)
	colnames(resmat) <- c("index", "loglik", "AIC", "conv", "upperbound", "lambda0", "lambda1", "mu0", "mu1", "q01", "q10", "pval_lrt")
	
	zf <- function(x){
		if (class(x) != "try-error"){
			return(TRUE)
		}else{
			return(FALSE)
		}
	}
	
	if (zf(fitx6) & zf(fitx5lv) & zf(fitx5mv) & zf(fitx4)){
		
		x <- fitx6
		resmat[1, 1:5] <- c(index, x$loglik, x$AIC, x$convergence, lpb)
		resmat[1, names(x$pars)] <- x$pars
		x <- fitx5lv
		resmat[2, 1:5] <- c(index, x$loglik, x$AIC, x$convergence, lpb)
		resmat[2, names(x$pars)] <- x$pars
		x <- fitx5mv
		resmat[3, 1:5] <- c(index, x$loglik, x$AIC, x$convergence, lpb)
		resmat[3, names(x$pars)] <- x$pars
		x <- fitx4
		resmat[4, 1:5] <- c(index, x$loglik, x$AIC, x$convergence, lpb)
		resmat[4, names(x$pars)] <- x$pars
		
		# lrt vs state-independent
		dpvec <- c(2,1,1)
		for (i in 1:3){
			ll <- 2 * (resmat[i,"loglik"] - resmat[4,"loglik"])
			if (ll < 0 ){
				resmat[i,"pval_lrt"] <- 1
			}else{
				resmat[i,"pval_lrt"] <- 1 - pchisq(ll, df=dpvec[i])
			}
				
		}		
		
	}else{
		resmat[1:4,1] <- index
		resmat[1:4, 2:12] <- NA		
	}
	

 
	return(resmat)	
}






DR_statistic <- function(x, return.mean = FALSE){
	
	rootnode <- length(x$tip.label) + 1
	
	sprates <- numeric(length(x$tip.label))
	for (i in 1:length(sprates)){
		node <- i
		index <- 1
		qx <- 0
		while (node != rootnode){
			el <- x$edge.length[x$edge[,2] == node]
			node <- x$edge[,1][x$edge[,2] == node]
			
			qx <- qx + el* (1 / 2^(index-1))
			
			index <- index + 1
		}
		sprates[i] <- 1/qx
	}
	
	if (return.mean){
		return(mean(sprates))		
	}else{
		names(sprates) <- x$tip.label
		return(sprates)
	}

}

# Simulate trees, conditional on getting at least minf frequency of 
#	the derived character state
simTreeBiSSE <- function(pars, tmax = 50, minf = 0, maxf = 1, nmin = 50, n_fixed=NULL){
	
	MAXBAD <- 500
 	badcount <- 0
 	
	while (1){
		
		if (is.null(n_fixed)){
			tmp <- tree.bisse(pars, max.t=tmax, x0=0)	
		}else{
			tmp <- tree.bisse(pars, max.taxa = n_fixed, x0=0)
			nmin <- n_fixed
		}
 
		if (!is.null(tmp)){
			ff <- sum(tmp$tip.state == 1) / length(tmp$tip.state)
			if (ff > minf & ff < maxf & length(tmp$tip.label) >= nmin){
				break
			}else if (badcount > MAXBAD){
				stop("Parameters cannot generate good tree/char combination\n")
			}
			badcount <- badcount + 1
		}
	}
	return(tmp)
}

asPhyloBaseClass <- function(phy){
	
	obj <- list(edge=phy$edge, Nnode = phy$Nnode, tip.label = phy$tip.label, edge.length = phy$edge.length)
	class(obj) <- "phylo"
	return(obj)
}


countParsimonyChanges <- function(phy, s){
 	
	phy <- asPhyloBaseClass(phy)
	
	if (length(unique(s)) == 1){
		return(0)
	}
	
	if (sum(names(s) == phy$tip.label) != length(phy$tip.label)){
		stop("Mismatch: state names vs tip labels\n")
	}
	
	# This fails occasionally.
	#phdata <- phyDat(s, type = "USER", levels = c(0,1))
 
	chars <- c("A", "C", "T", "G")
	if (sum(s  == 0) > 0){
		s <- s + 1
	}
	
	cvec <- as.matrix(chars[s], ncol=1) 
	rownames(cvec) <- names(s)
	phdata <- phyDat(cvec, type = "DNA")

	return(parsimony(phy, phdata ))
 
}

getWeightedDR <- function(phy, s, DRvec){
 	
	phy <- asPhyloBaseClass(phy)
	
	if (length(unique(s)) == 1){
		return(0)
	}
	
	if (sum(names(s) == phy$tip.label) != length(phy$tip.label)){
		stop("Mismatch: state names vs tip labels\n")
	}
	
	# This fails occasionally.
	#phdata <- phyDat(s, type = "USER", levels = c(0,1))
 
	chars <- c("A", "C", "T", "G")
	if (sum(s  == 0) > 0){
		s <- s + 1
	}
	
	cvec <- as.matrix(chars[s], ncol=1) 
	rownames(cvec) <- names(s)
	phdata <- phyDat(cvec, type = "DNA")

	at <- acctran(phy, phdata)

	rootnode <- length(phy$tip.label) + 1

	curr_state <- s[phy$tip.label[1]]
	which_edge <- which(phy$edge[,2] == 1)
	parent <- phy$edge[which_edge, 1]
	while(parent != rootnode){
		if (at$edge.length[which_edge] == 1){
			curr_state <- abs(curr_state - 1)
		}
		which_edge <- which(phy$edge[,2] == parent)
		parent <- phy$edge[which_edge, 1]
		
	}
	at$rootstate <- curr_state
	
	tipset <- phy$tip.label
	
	ratesums      <- numeric(2)
	countsums     <- numeric(2)
	
	# Do this directly w vcv matrix:
	vv <- cophenetic.phylo(at)
	
	# First get rid of all terminal changes:
	
	while (length(tipset) > 0){
		curtip <- which(rownames(vv) == tipset[1])
		in_group <- which(vv[,curtip] == 0)
		
		tips <- rownames(vv)[in_group]
		
		index <- s[tipset[1]]
		
		ratesums[ index  ] <- ratesums[ index] + sum(DRvec[tips]) / length(tips)
		countsums[ index] <- countsums[ index  ] + 1
		tipset <- setdiff(tipset, tips)	
  
	}
	
 
	ll <- list(rates = ratesums/countsums, ratesums = ratesums, countsums = countsums)
 	return(ll)
 	
}

 

simTraitsParsimonyCriterion <- function(phy, par, changes, tol, fail_tol = 5000){
 
	states <- 0
	
	changes_min <- changes - round(changes * tol)
	changes_max <- changes + round(changes * tol)
	
	counter <- 0
	
	while(1){
		
		states <- sim.character(phy, par, x0=sample(c(1,0), 1), model="mk2")
		cx <- countParsimonyChanges(phy, states)
		if (cx >= changes_min & cx <= changes_max){
			break
		}
		counter <- counter + 1
		
		if (counter > fail_tol){
			cat("Too many failed attempts in simTraitsParsimonyCriterion\n")
			cat("It appears that transition rates are mismatched to tree and trait data\n")
			stop("See tol and fail_tol arguments to this function")
		}
		
	}
 
	return(states)	
}

fitMK1 <- function(phy, states, mx = 10)
{
	lfx <- make.mk2(phy, states)
	lfx <- constrain(lfx, formulae = list(q01 ~ q10))
	qres <- optimize(lfx, interval=c(0.00001, mx), maximum=T)
	return(qres$maximum)
}



fitMK2 <- function(phy, states, nopt = 1){
	
	lfx <- make.mk2(phy, states)
	lfx2 <- function(pars){
		pars <- exp(pars)
		return( lfx(pars) )
	}
	
	for (i in 1:nopt){
		
		badcount <- 0
		iv <- log(runif(2, 0, 0.5))
		resx <- try(optim(iv , lfx2, control=list(fnscale=-1)))
		while (class(resx) == 'try-error'){
			iv <- log(runif(2, 0, 0.5))
			resx <- try(optim(iv , lfx2, control=list(fnscale=-1)))
			
			badcount <- badcount + 1
			if (badcount > MAXBAD){
				stop("Too many fails in fitDiversitree\n")
			}
		}
		
		if (i == 1){
			best <- resx
		}else{
			if (best$value < resx$value){
				best <- resx
			}
		}
		
	}
	
	
	fres <- list(pars=exp(best$par), loglik=best$value)
	fres$AIC <- -2*fres$loglik + 2*2
	fres$counts <- best$counts
 
	fres$convergence <- best$convergence
	fres$message <- best$message
	return(fres)
	
}

# Fast Intuitive analysis of State-dependent Speciation Extinction rates
#      incomplete = TRUE  : allows trait data to be subset of tree, e.g., incomplete sampling for trait data
#                  
FISSE.binary <- function(phy, states, reps = 1000, tol=0.1, qratetype = "mk", incomplete = TRUE, ...){
 
	mism <- setdiff(names(states), phy$tip.label)
	if (length(mism) > 0){
		stop("Error: Trait data includes taxa that are not present in tree\n")
	}
	
	if (! incomplete){
		if (length(intersect(names(states), phy$tip.label)) != length(phy$tip.label) ){
			stop("error in names matching between tree tips and state vector")
		}		
	}
	
	# The ES measures:
	dx <- DR_statistic(phy)	
	
	# now drop tree to same set of taxa in states dataset:
	phy <- drop.tip(phy, tip = setdiff(phy$tip.label, names(states)))
 	
 	dx <- dx[phy$tip.label]
 
	states <- states[phy$tip.label] 

	lam0 <- dx[names(states)[states == 0] ]
	lam1 <- dx[names(states)[states == 1] ]
 
	nc <-  countParsimonyChanges(phy, states)	
	
	qq <- 0
	
	if (qratetype == "parsimony"){
		qq <- nc / sum(phy$edge.length)		
	}else if (qratetype == "mk"){
		qq <- fitMK1(phy, states)
	}else{
		stop("Unsupported or invalid option to estimate q")
	}
	
	
	#mk2res <- fitMK2(phy, states)

	
	mm <- matrix(NA, nrow=reps, ncol=2)
	for (i in 1:reps){
	
		#tset <- simTraitsParsimonyCriterion(phy, mk2res$pars, nc)
		tset <- try(simTraitsParsimonyCriterion(phy, c(qq,qq), nc, tol, ...))
		if (class(tset) != "try-error"){
			mm[i,1] <- mean(dx[names(tset)[tset == 0]])
			mm[i,2] <- mean(dx[names(tset)[tset == 1]])			
		}else{
			obj <- list(lambda0 = mean(lam0), lambda1 = mean(lam1), pval = NA, null_mean_diff = NA, null_sd = NA, nchanges_parsimony = nc, qpars = qq)
			return(obj)	
		}
 
	}
	
	delta <- mm[,2] - mm[,1]
	delta_true <- mean(lam1) - mean(lam0)
	pval <- sum(delta_true > delta) / (reps + 1)
	
	obj <- list(lambda0 = mean(lam0), lambda1 = mean(lam1), pval = pval, null_mean_diff = mean(delta), null_sd = sd(delta), nchanges_parsimony = nc, qpars = qq)
	return(obj)
}


simulateCharacter <- function(phy, qval, minf){
	
 	maxbadcount <- 0
	good <- FALSE
	while (!good){
		if (maxbadcount > 25000){
			stop("maxbadcount exceeded\n")
		}
		
		states <- sim.character(phy, c(qval, qval), x0 = sample(c(0,1), 1), model = "mk2")
 
		states <- states[phy$tip.label]
		tx <- table(states)
		if (length(tx) > 1 & (min(tx) > minf*length(phy$tip.label))){
			good <- TRUE
		}
 			
		maxbadcount <- maxbadcount + 1
	}
 
	return(states)
} 
 

########################## 
# Batch processing function

runAnalyses <- function(treefile, traitfile, id = "xxx", bisse_opt = 5){
	
	v <- read.tree(treefile)
	t1 <- read.csv(traitfile, header=F, stringsAsFactors=F)
	states <- t1[,2]
	names(states) <- t1[,1]
	
	ff <- FISSE.binary(v, states, reps=2000)
	bisse <- fitDiversitree_allmodels(v, states, nopt=bisse_opt)
	
	ntips = length(v$tip.label)
	f0 <- sum(states == 0) / ntips
	f1 <- sum(states == 1) / ntips
	
	
	
	fisse <- c(ntips=ntips, f0=f0, f1=f1, unlist(ff))
	res <- list(id=id, traits=traitfile, tree=treefile, fisse=fisse, bisse=bisse)
	return(res)
}



# Fast Intuitive analysis of State-dependent Speciation Extinction rates
#  Weighted version
#  Rates by state are weighted such that each parsimony origin of a state counts as a single 
#    point, regardless of how many tips are included.
#
#  This function was not used in Rabosky & Goldberg, Evolution, 2017
#  Has relatively low power but is more robust to phylogenetic pseudoreplication
#   
FISSE.weighted.binary <- function(phy, states, reps = 1000, tol=0.1, qratetype = "mk", ...){
	
	if (length(intersect(names(states), phy$tip.label)) != length(phy$tip.label) ){
		stop("error in names matching between tree tips and state vector")
	}
	states <- states[phy$tip.label] 
	
	# The ES measures:
	dx <- DR_statistic(phy)
	tmp <- getWeightedDR(phy, states, dx)
	lam0 <- tmp$rates[1]
	lam1 <- tmp$rates[2]
 
	nc <-  countParsimonyChanges(phy, states)	
	
	qq <- 0
	
	if (qratetype == "parsimony"){
		qq <- nc / sum(phy$edge.length)		
	}else if (qratetype == "mk"){
		qq <- fitMK1(phy, states)
	}else{
		stop("Unsupported or invalid option to estimate q")
	}
	
	
	#mk2res <- fitMK2(phy, states)

	
	mm <- matrix(NA, nrow=reps, ncol=2)
	for (i in 1:reps){
	
		#tset <- simTraitsParsimonyCriterion(phy, mk2res$pars, nc)
		tset <- try(simTraitsParsimonyCriterion(phy, c(qq,qq), nc, ...))
		if (class(tset) != "try-error"){
			
			tmp <- getWeightedDR(phy, tset, dx)
			
			mm[i,1] <- tmp$rates[1]
			mm[i,2] <- tmp$rates[2]		
		}else{
			obj <- list(lambda0 = lam0, lambda1 = lam1, pval = NA, null_mean_diff = NA, null_sd = NA, nchanges_parsimony = nc, qpars = qq)
			return(obj)	
		}
 
	}
	
	delta <- mm[,2] - mm[,1]
	delta_true <- lam1 - lam0
	pval <- sum(delta_true > delta) / (reps + 1)
	
	obj <- list(lambda0 = lam0, lambda1 = lam1, pval = pval, null_mean_diff = mean(delta), null_sd = sd(delta), nchanges_parsimony = nc, qpars = qq)
	return(obj)
}














