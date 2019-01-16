# R software implementing Patient-Net (P-Net): 
# Network-based ranking of patients with respect to a given phenotype/outcome. 

##*********************************************************************##
# Copyright (C) 2017 Jessica Gliozzo

# This file is part of P-Net library. 

# P-Net is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
##*********************************************************************##

library(RANKS);
library(PerfMeas);


###################################################################################
# Methods implementing:
# P-net algorithm implemented through loo, cv and heldout. 
# Automatic setting of net threshold is implemented through efficient internal loo.
###################################################################################

# Function to select a set of thresholds for a network. 
# Thresholds corresponding to different quantiles are selected. Minima of x are deleted to avoid distributions skewed towards 0 (the expected minimum).
# Input:
# x : vector or matrix of the values.
# probs : the probabilities to compute the corresponding quantiles.
# Output:
# a vector with the computed thresholds.
compute.thresholds <- function (x, probs) {
  x <- x[x>min(x)]; # minima are not considered
  thresh <- quantile(x, probs=probs);
  names(thresh)<-probs;
  return(thresh);
}


# Method to find the "optimal" filtering threshold by loo.
# The loo is used to select the optimal network threshold:
# The optimal threshold in the loo is selected using the AUC as the objective function.
# Input:
# W : a symmetric matrix representing similarity between examples. It can be a kernel or a correlation matrix, or any other matrix representing meaningful similarities between examples.
# ind.pos : integer vector with the indices of the positive examples. Indices corresponding to the row (column) of W.
# ind.test : integer vector with the indices of the examples to be tested (def: 1:nrow(W)). Note the default corresponds to the classical loo procedure.
# score : Local score function to be considered. It may be one of the following:
#                 - eav.score (default)
#                 - NN.score
#                 - KNN.score 
#                 - WSLD score
#                 - tot.score
# 		  		  - diff.score
#		          - dnorm.score
# probs: probabilities relative to the quantiles to be computed (def. 0.1 to 1 by 0.1 steps). Minimum values of W (usually zeros) are not considered.
#        Note: the numeric vector must be monotonically crescent.
# ... : optional arguments for the function score.
# Output:
# a list with the following elements:
#   - quantile : the prob of the selected quantile.
#   - thresh : the threshold corresponding to the selected quantile.
#   - auc : the AUC computed using as threshold the "optimal" threshold thresh.
setGeneric("optimize.thresh.by.loo", 
                 function(W, ind.pos, ind.test=c(1:nrow(W)), score=eav.score, probs=seq(0.1, 1, 0.1), ...) standardGeneric("optimize.thresh.by.loo"));
				 
setMethod("optimize.thresh.by.loo", signature(W="matrix"),
  function(W, ind.pos, ind.test=c(1:nrow(W)), score=eav.score, probs=seq(0.1, 1, 0.1), ...) {
    diag(W) <- 0; # this is for the loo in 1 run of the algorithm
    opt.auc <- 0;
    opt.q <- 0;
    opt.thresh <- 0;
	n <- nrow(W);
	targets <- integer(n);  
    targets[ind.pos] <- 1;
	
	thresh <- compute.thresholds(W, probs);
	n.thresh <- length(thresh);
    W.filtered <- W;
    for (i in 1:n.thresh) {
	  W.filtered[W.filtered<thresh[i]] <- 0;
	  scores <- score(W.filtered, ind.test, ind.pos, ...);
	  res <- AUC.single(scores, targets[ind.test]);
	  if (res >= opt.auc) {
	    opt.auc <- res;
		opt.q<-as.numeric(names(thresh)[i]);
		opt.thresh<-thresh[i];	  
	  }
	}  
  
    return(list(quantile=opt.q, thresh=opt.thresh, auc=opt.auc)); 
})

# Method that implements the P-Net algorithm using a double loo.
# An external loo is performed to evaluate the performance of the method, and an internal loo is performed to select the optimal threshold.
# The threshold to filter the net for the external loo is obtained by averaging across the thresholds estimated at each internal loo.
# Note that if intloo=FALSE no internal loo is executed and the net is filtered according to the quantile of the argument probs.
# Inputs:
# W : a symmetric matrix representing similarity between examples. It can be a kernel or a correlation matrix, or any other matrix representing meaningful similarities between examples.
# ind.pos : integer vector with the indices of the positive examples. Indices corresponding to the row (column) of W.
# score : Local score function to be considered. It may be one of the following:
#                 - eav.score (default)
#                 - NN.score
#                 - KNN.score 
#                 - WSLD score
#                 - tot.score
# 		  		  - diff.score
#		          - dnorm.score
# probs: probabilities relative to the quantiles to be computed (def. 0.1 to 1 by 0.1 steps). Minimum values of W (usually zeros) are not considered.
#        Note: the numeric vector must be monotonically crescent. If intloo=FALSE, probs represents the quantile that will be used to filter the network.
# intloo : logical. If TRUE (default) an internal loo is performed to find the optimal threshold, otherwise no internal loo is performed and the probs value is used to compute the quantile for filtering.
# ... : optional arguments for the function score.
# Output:
# a list with the following elements:
#   - s : a vector with scores computed by loo.
#   - q : the quantile vector corresponding to the quantiles selected for each patient.
setGeneric("pnet.loo", 
                 function(W, ind.pos, score=eav.score, probs=seq(0.1, 1, 0.1), intloo=TRUE, ...) standardGeneric("pnet.loo"));
				 
setMethod("pnet.loo", signature(W="matrix"),
  function(W, ind.pos, score=eav.score, probs=seq(0.1, 1, 0.1),  intloo=TRUE, ...) {
    m <- nrow(W);
	s  <- numeric(m);
	names(s) <- rownames(W); 
	if (intloo) {
	  q <-numeric(m);
	  names(q) <- rownames(W); 
	  # computing the optimal threshold for each row of W (by internal loo)
	  for (i in 1:m)  {
	    ind.test <- setdiff(1:m,i);
	    q[i] <- optimize.thresh.by.loo(W, setdiff(ind.pos,i), ind.test=ind.test, score=score, probs=probs, ...)$quantile;
	  }
	  thresh <- quantile(W,probs=mean(q));
	} else {
	  thresh <- quantile(W,probs=probs);
	  q=probs;
	}
	
	W[W<thresh] <- 0;
	diag(W) <- 0;  # efficient external loo
	s <- score(W, 1:m, ind.pos, ...);
	mm<-max(s);
	if (mm > 0)
	  s <- s / mm;
    return (list(s=s, q=q)); 
})


# Method that implements the P-Net algorithm using a double loo (selection of an optimal threshold for each patient).
# An external loo is performed to evaluate the performance of the method, and an internal loo is performed to select the optimal threshold
# for each patient. A different threshold is selected for each patient.
# Note that if intloo=FALSE no internal loo is executed and the net is filtered according to the quantile of the argument probs.
# Inputs:
# W : a symmetric matrix representing similarity between examples. It can be a kernel or a correlation matrix, or any other matrix representing meaningful similarities between examples.
# ind.pos : integer vector with the indices of the positive examples. Indices corresponding to the row (column) of W.
# score : Local score function to be considered. It may be one of the following:
#                 - eav.score (default)
#                 - NN.score
#                 - KNN.score 
#                 - WSLD score
#                 - tot.score
# 		          - diff.score
#		          - dnorm.score
# probs: probabilities relative to the quantiles to be computed (def. 0.1 to 1 by 0.1 steps). Minimum values of W (usually zeros) are not considered.
#        Note: the numeric vector must be monotonically crescent. If intloo=FALSE, probs represents the quantile that will be used to filter the network.
# intloo : logical. If TRUE (default) an internal loo is performed to find the optimal threshold, otherwise no internal loo is performed and the probs value is used to compute the quantile for filtering.
# ... : optional arguments for the function score.
# Output:
# a list with the following elements:
#   - s : a vector with scores computed by loo.
#   - q : the quantile vector corresponding to the quantiles selected for each patient.
setGeneric("pnet.loo.per.patient.filtering", 
                 function(W, ind.pos, score=eav.score, probs=seq(0.1, 1, 0.1), intloo=TRUE, ...) standardGeneric("pnet.loo.per.patient.filtering"));
				 
setMethod("pnet.loo.per.patient.filtering", signature(W="matrix"),
  function(W, ind.pos, score=eav.score, probs=seq(0.1, 1, 0.1), intloo=TRUE, ...) {
    m <- nrow(W);
	s  <- numeric(m);
	names(s) <- rownames(W); 
	if (intloo) {
	  q <-numeric(m);
	  names(q) <- rownames(W); 
	# computing the optimal threshold for each row of W (by internal loo), and filtering the corresponding row with the computed optimal threshold
	  for (i in 1:m)  {
	    ind.test <- setdiff(1:m,i);
	    q[i] <- optimize.thresh.by.loo(W, setdiff(ind.pos,i), ind.test=ind.test, score=score, probs=probs, ...)$quantile;
	    thresh <- quantile(W[i,],probs=q[i]);	
	    W[i,W[i,]<thresh] <- 0; # by row filtering
	  }
	} else {
	  thresh <- quantile(W,probs=probs);
	  q=probs;
	  W[W<thresh] <- 0;
	}
	diag(W) <- 0;  # efficient external loo
	s <- score(W, 1:m, ind.pos, ...);
    return (list(s=s, q=q)); 
})


# Method that implements the P-Net algorithm using cross-validation.
# An external CV is performed to evaluate the performance of the method, and an internal loo is performed to select the optimal threshold.
# The threshold to filter the net for the external cv is obtained by averaging across the thresholds estimated at each internal loo
# for each fold.
# Note that if intloo=FALSE no internal loo is executed and the net is filtered according to the quantile of the argument probs.
# W : a symmetric matrix representing similarity between examples. It can be a kernel or a correlation matrix, or any other matrix representing meaningful similarities between examples.
# ind.pos : integer vector with the indices of the positive examples. Indices corresponding to the row (column) of W.
# kk : number of folds (def. 5)
# score : Local score function to be considered. It may be one of the following:
#                 - eav.score (default)
#                 - NN.score
#                 - KNN.score 
#                 - WSLD score
#                 - tot.score
# 		          - diff.score
#		          - dnorm.score
# probs: probabilities relative to the quantiles to be computed (def. 0.1 to 1 by 0.1 steps). Minimum values of W (usually zeros) are not considered.
#        Note: the numeric vector must be monotonically crescent. If intloo=FALSE, probs represents the quantile that will be used to filter the network.
# intloo : logical. If TRUE (default) an internal loo is performed to find the optimal threshold, otherwise no internal loo is performed and the probs value is used to compute the quantile for filtering.
# seed : init seed for the random construction of the folds (def. 1).
# ... : optional arguments for the function score.
# Output:
# a list with the following elements:
#   - s : a vector with scores computed by k-fold cross validation.
#   - q : the quantile vector corresponding to the quantiles selected at each fold.
setGeneric("pnet.cv", 
                 function(W, ind.pos, kk=5, score=eav.score, probs=seq(0.1, 1, 0.1), intloo=TRUE, seed=NULL, ...) standardGeneric("pnet.cv"));
				 
setMethod("pnet.cv", signature(W="matrix"), 
  function(W, ind.pos, kk=5, score=eav.score, probs=seq(0.1, 1, 0.1), intloo=TRUE, seed=NULL, ...) {
	score.name <- return.name(score);
	m <- nrow(W);
	x <- 1:m;
	f <- do.stratified.cv.data(1:m, ind.pos, k=kk, seed=seed);    
	s  <- numeric(m);
	names(s) <- rownames(W);
	if (intloo) { 
	  q <- numeric(kk);
	  # computing the optimal threshold for each fold of W (by internal loo)
	  for (i in 1:kk)  {
	    ind.train <- setdiff(x, c(f$fold.non.positives[[i]],f$fold.positives[[i]]));
	    ind.pos.train<-setdiff(ind.pos,f$fold.positives[[i]]);
	    if (score.name == "KNN.score")
	            q[i] <- optimize.thresh.by.loo(W, ind.pos.train, ind.train, score=score, probs=probs, norm=FALSE, k=k, ...)$quantile
	    else if (score.name == "WSLD.score")								   
		    q[i] <- optimize.thresh.by.loo(W, ind.pos.train, ind.train, score=score, probs=probs, norm=FALSE, d=d, ...)$quantile
	    else if (score.name == "tot.score" || score.name == "diff.score" || score.name == "dnorm.score")			   
		    q[i] <- optimize.thresh.by.loo(W, ind.pos.train, ind.train, score=score, probs=probs)$quantile
	    else
		    q[i] <- optimize.thresh.by.loo(W, ind.pos.train, ind.train, score=score, probs=probs, norm=FALSE, ...)$quantile;

	  }
	  thresh <- quantile(W,probs=mean(q));	
	} else {
	  thresh <- quantile(W,probs=probs);
	  q=probs;	
	}
	W[W<thresh] <- 0;
	# external cross validation
	for (i in 1:kk)  {
	  ind.test <- c(f$fold.non.positives[[i]],f$fold.positives[[i]]);
	  ind.pos.fold<-setdiff(ind.pos,f$fold.positives[[i]]);
	  if (score.name == "KNN.score")  				  
 	  	s[ind.test] <- score(W, ind.test, ind.pos.fold, k=k, norm=FALSE, ...)   
    	  else if (score.name == "WSLD.score")		  
 	  	s[ind.test] <- score(W, ind.test, ind.pos.fold, d=d, norm=FALSE, ...) 
    	  else if (score.name == "tot.score" || score.name == "diff.score" || score.name == "dnorm.score") 
          	s[ind.test] <- score(W, ind.test, ind.pos.fold) 
    	  else										  
 	  	s[ind.test] <- score(W, ind.test, ind.pos.fold, norm=FALSE, ...); 
	}
	mm<-max(s);
	if (mm > 0)
	  s <- s / mm;
	return (list (s=s, q=q)); 
})


# Method that implements the P-Net algorithm using held-out.
# A loo on the training set is performed to select the optimal threshold.
# The threshold to filter the net for the external held-out is obtained by averaging across the thresholds estimated at each internal loo
# for each fold.
# Note that if intloo=FALSE no internal loo is executed and the net is filtered according to the quantile of the argument probs.
# Input:
# W : a symmetric matrix representing similarity between examples. It can be a kernel or a correlation matrix, or any other matrix representing meaningful similarities between examples.
# ind.pos : integer vector with the indices of the positive examples. Indices correspond to the row (column) of W.
# test : ratio of the test data (def. 0.3). I.e. if test is 0.3 of the available data, the train set will be the remaining 0.7.
#         Alternatively may be a vector with the indices of the test data. Indices corrrespond to the rows of W.
# score : Local score function to be considered. It may be one of the following:
#                 - eav.score (default)
#                 - NN.score
#                 - KNN.score 
#                 - WSLD score
#                 - tot.score
# 		          - diff.score
#		          - dnorm.score
# probs: probabilities relative to the quantiles to be computed (def. 0.1 to 1 by 0.1 steps). Minimum values of W (usually zeros) are not considered.
#        Note: the numeric vector must be monotonically crescent. If intloo=FALSE, probs represents the quantile that will be used to filter the network.
# intloo : logical. If TRUE (default) an internal loo is performed to find the optimal threshold, otherwise no internal loo is performed and the probs value is used to compute the quantile for filtering.
# seed : init seed for the random construction of the folds (def. 1).
# ... : optional arguments for the function score.
# Output:
# a list with the following elements:
#   - s : a vector of the scores (including both training and test data).
#   - test : indices of the elements of the test set (indices of the vector s).
#   - q : the quantile  selected on the training set.
setGeneric("pnet.heldout", 
                 function(W, ind.pos, test=0.3, score=eav.score, probs=seq(0.1, 1, 0.1), intloo=TRUE, seed=NULL, ...) standardGeneric("pnet.heldout"));
				 
setMethod("pnet.heldout", signature(W="matrix"),
  function(W, ind.pos, test=0.3, score=eav.score, probs=seq(0.1, 1, 0.1), intloo=TRUE, seed=NULL, ...) {
    set.seed(seed);
	score.name <- return.name(score);
	m <- nrow(W);
	x <- 1:m;
	s  <- numeric(m);
	names(s) <- rownames(W); 
	if (is.integer(test[1]))
	   ind.test <- test
	else
	   ind.test <- sample(x, trunc(m*test));   
	ind.train <- setdiff(x, ind.test);
	ind.pos.train <- setdiff(ind.pos,ind.test); 
	
	# computing the optimal threshold by loo w.r.t. the training set	
	if (intloo) {
	  if (score.name == "KNN.score") 
	     q <- optimize.thresh.by.loo(W, ind.pos.train, ind.train, score=score, probs=probs, norm=FALSE, k=k, ...)$quantile
	  else if (score.name == "WSLD.score")								   
             q <- optimize.thresh.by.loo(W, ind.pos.train, ind.train, score=score, probs=probs, norm=FALSE, d=d, ...)$quantile
	  else if (score.name == "tot.score" || score.name == "diff.score" || score.name == "dnorm.score")			   
	     q <- optimize.thresh.by.loo(W, ind.pos.train, ind.train, score=score, probs=probs)$quantile
	  else
	     q <- optimize.thresh.by.loo(W, ind.pos.train, ind.train, score=score, probs=probs, norm=FALSE, ...)$quantile;
          }
	else {
	   q <- probs;
        }
	thresh <- quantile(W,probs=q);	
	W[W<thresh] <- 0;
	# test on held.out data
        if (score.name == "KNN.score")  				  
 	   s <- score(W, x, ind.pos.train, k=k, norm=FALSE, ...)   
    	else if (score.name == "WSLD.score")		  
 	   s <- score(W, x, ind.pos.train, d=d, norm=FALSE, ...) 
    	else if (score.name == "tot.score" || score.name == "diff.score" || score.name == "dnorm.score") 
           s <- score(W, x, ind.pos.train) 
    	else										  
 	   s <- score(W, x, ind.pos.train, norm=FALSE, ...);

	mm<-max(s);
	if (mm > 0)
	  s <- s / mm;
    return (list(s=s, test=ind.test, q=q)); 
})


#######################################################################
# Methods implementing:
# Non parametric test to validate patients' ranking and biomarkers sets
#######################################################################

##################################################
# Statistical tests for patients' ranking
##################################################


# Function that implements a non parametric test to assess the statistical significance of the overall ranking obtained with P-Net through loo.
# The test is based on multiple permutations of the labels of the examples, and a p.value is estimated by comparing the AUC 
# achieved with the true labels with those obtained by random permutations.
# The threshold to filter the net for the external loo is obtained by averaging across the thresholds estimated at each internal loo.
# Note that if intloo=FALSE no internal loo is executed and the net is filtered according to the quantile of the argument probs.
# Input:
# W : a symmetric matrix representing similarity between examples. It can be a kernel or a correlation matrix, or any other matrix representing meaningful similarities between examples. 
# ind.pos : integer vector with the indices of the positive examples. Indices correspond to the row (column) of W.
# score : Local score function to be considered. It may be one of the following:
#                 - eav.score (default)
#                 - NN.score
#                 - KNN.score 
#                 - WSLD score
# 		          - tot.score
# 		          - diff.score
#		          - dnorm.score
# probs: probabilities relative to the quantiles to be computed (def. 0.1 to 1 by 0.1 steps). Minimum values of W (usually zeros) are not considered.
#        Note: the numeric vector must be monotonically crescent. If intloo=FALSE, probs represents the quantile that will be used to filter the network.
# rep : number of the repetitions of the label permutations (def: 1000).
# intloo : logical. If TRUE (default) an internal loo is performed to find the optimal threshold, otherwise no internal loo is performed and the probs value is used to compute the quantile for filtering.
# seed : init seed for the random permutations (def. NULL), If NULL no initialization is performed.
# ... : optional arguments for the function score.
# Output:
# a list with the following elements:
#   - rank.p.value : the computed p-value for the overall ranking.
#   - auc : value of the computed AUC by loo.
#   - scores : vector of the scores computed for each gene.
ranking.assessment.loo <-function(W, ind.pos, score=eav.score, probs=seq(0.1, 1, 0.1), rep=1000, intloo=TRUE, seed=NULL, ...) {

 set.seed(seed);
 m <- nrow(W);
 # pnet computed on the true labels
 res <- pnet.loo(W, ind.pos, score=score, probs=probs,  intloo=intloo, ...);
 s <- res$s;
 names(s) <- rownames(W); 
 q <- mean(res$q);
 target <- numeric(m);
 target[ind.pos] <- 1;
 auc <- AUC.single(s, target);
 
 # random permutations
 n.pos<-length(ind.pos);
 x <- 1:m;
 rank.p.value <- 0.0;
 for (i in 1:rep) {
   ind.pos.random <- sample(x, n.pos);
   target <- rep(0,m);
   target[ind.pos.random] <- 1;
   res <- pnet.loo(W, ind.pos.random, score=score, probs=q,  intloo=FALSE, ...);
   auc.random <- AUC.single(res$s, target);
   if (auc.random >= auc)
     rank.p.value <- rank.p.value + 1;
 }
 return(list(rank.p.value=rank.p.value/rep, auc=auc, scores=s));
}

# Function that implements a non parametric test to assess the statistical significance of the overall ranking obtained with P-Net through cross-validation.
# The test is based on multiple permutations of the labels of the examples, and a p.value is estimated by comparing the AUC 
# achieved with the true labels with those obtained by random permutations.
# The threshold to filter the net for the external cross-validation is obtained by averaging across the thresholds estimated at each internal loo.
# Note that if intloo=FALSE no internal loo is executed and the net is filtered according to the quantile of the argument probs.
# Input:
# W : a symmetric matrix representing similarity between examples. It can be a kernel or a correlation matrix, or any other matrix representing meaningful similarities between examples.
# ind.pos : integer vector with the indices of the positive examples. Indices correspond to the row (column) of W.
# kk : number of folds (def. 5).
# score : Local score function to be considered. It may be one of the following:
#                 - eav.score (default)
#                 - NN.score
#                 - KNN.score 
#                 - WSLD score
# 		          - tot.score
# 		          - diff.score
#		          - dnorm.score
# probs: probabilities relative to the quantiles to be computed (def. 0.1 to 1 by 0.1 steps). Minimum values of W (usually zeros) are not considered.
#        Note: the numeric vector must be monotonically crescent. If intloo=FALSE, probs represents the quantile that will be used to filter the network.
# rep : number of the repetitions of the label permutations (def: 1000).
# intloo : logical. If TRUE (default) an internal loo is performed to find the optimal threshold, otherwise no internal loo is performed and the probs value is used to compute the quantile for filtering.
# seed : init seed for the random permutations (def. NULL). If NULL no initialization is performed.
# ... : optional arguments for the function score.
# Output:
# a list with the following elements:
#   - rank.p.value : the computed p-value for the overall ranking.
#   - auc : value of the computed AUC by cross-validation.
#   - scores : vector of the scores computed for each gene.
ranking.assessment.cv <-function(W, ind.pos, kk=5, score=eav.score, probs=seq(0.1, 1, 0.1), rep=1000, intloo=TRUE, seed=NULL, ...) {

 set.seed(seed);
 m <- nrow(W);
 # pnet computed on the true labels
 res <- pnet.cv(W, kk=kk, ind.pos=ind.pos, score=score, probs=probs, intloo=intloo, seed=seed, ...)
 s <- res$s;
 names(s) <- rownames(W); 
 q <- mean(res$q);
 target <- numeric(m);
 target[ind.pos] <- 1;
 auc <- AUC.single(s, target);
 
 # random permutations
 n.pos<-length(ind.pos);
 x <- 1:m;
 rank.p.value <- 0.0;
 for (i in 1:rep) {
   ind.pos.random <- sample(x, n.pos);
   target <- rep(0,m);
   target[ind.pos.random] <- 1;
   res <- pnet.cv(W, kk=kk, ind.pos=ind.pos.random, score=score, probs=q,  intloo=FALSE, ...);
   auc.random <- AUC.single(res$s, target);
   if (auc.random >= auc)
     rank.p.value <- rank.p.value + 1;
 }
 return(list(rank.p.value=rank.p.value/rep, auc=auc, scores=s));
}

# Function that implements a non parametric test to assess the statistical significance of the overall ranking obtained with P-Net through held out.
# The test is based on multiple permutations of the labels of the examples, and a p.value is estimated by comparing the AUC 
# achieved with the true labels with those obtained by random permutations.
# Note that if intloo=FALSE no internal loo is executed and the net is filtered according to the quantile of the argument probs.
# Input:
# W : a symmetric matrix representing similarity between examples. It can be a kernel or a correlation matrix, or any other matrix representing meaningful similarities between examples. 
# ind.pos : integer vector with the indices of the positive examples. Indices correspond to the row (column) of W.
# test : ratio of the test data (def. 0.3). I.e. if test is 0.3 of the available data, the train set will be the remaining 0.7.
#         Alternatively may be a vector with the indices of the test data. Indices correspond to the rows of W.
# score : Local score function to be considered. It may be one of the following:
#                 - eav.score (default)
#                 - NN.score
#                 - KNN.score 
#                 - WSLD score
# 		          - tot.score
# 		          - diff.score
#		          - dnorm.score
# probs: probabilities relative to the quantiles to be computed (def. 0.1 to 1 by 0.1 steps). Minimum values of W (usually zeros) are not considered.
#        Note: the numeric vector must be monotonically crescent. If intloo=FALSE, probs represents the quantile that will be used to filter the network.
# rep : number of the repetitions of the label permutations (def: 1000).
# intloo : logical. If TRUE (default) an internal loo is performed to find the optimal threshold, otherwise no internal loo is performed and the probs value is used to compute the quantile for filtering. 
# seed : init seed for the random permutations (def. NULL). If NULL no initialization is performed.
# ... : optional arguments for the function score. 
# Output:
# a list with the following elements:
#   - rank.p.value : the computed p-value for the overall ranking.
#   - auc : value of the computed AUC on the test data.
#   - scores : a vector of the scores (including both training and test data).
#   - test : indices of the elements of the test set (indices of the vector scores). 
ranking.assessment.heldout <-function(W, ind.pos, test=0.3, score=eav.score, probs=seq(0.1, 1, 0.1), rep=1000, intloo=TRUE, seed=NULL, ...) {

 set.seed(seed);
 m <- nrow(W);
 # pnet computed on the true labels
 res <- pnet.heldout(W, ind.pos=ind.pos, test=test, score=score, probs=probs, intloo=intloo, seed=seed, ...);
 s <- res$s;
 names(s) <- rownames(W); 
 q <- res$q;
 target <- numeric(m);
 target[ind.pos] <- 1;
 ind.test <- res$test;
 s.test<-s[ind.test];
 target.ind.test<-target[ind.test];
 auc <- AUC.single(s.test, target.ind.test);
 
 # random permutations
 n.pos<-length(ind.pos);
 x <- 1:m;
 rank.p.value <- 0.0;
 for (i in 1:rep) {
   ind.pos.random <- sample(x, n.pos);
   target <- rep(0,m);
   target[ind.pos.random] <- 1;
   res <- pnet.heldout(W, ind.pos=ind.pos.random, test=ind.test, score=score, probs=q, intloo=FALSE,  ...);
   s.test<-res$s[ind.test];
   target.ind.test<-target[ind.test];
   auc.random <- AUC.single(s.test, target.ind.test);
   if (auc.random >= auc)
     rank.p.value <- rank.p.value + 1;
 }
 return(list(rank.p.value=rank.p.value/rep, auc=auc, scores=s, test=ind.test));
}


##################################################
# Statistical tests for Biomarker set assessment
##################################################


# Function that implements a non parametric test to assess the statistical significance of a given biomarker set through loo.
# The test is based on multiple random resamplings of the biomarker set, and a p.value is estimated by comparing the AUC 
# achieved with the assessed biomarker set with those obtained by random resampling.
# The threshold to filter the net for the external loo is obtained by averaging across the thresholds estimated at each internal loo.
# Note that if intloo=FALSE no internal loo is executed and the net is filtered according to the quantile of the argument probs.
# Input:
# C : a matrix representing examples and the associated features. Rows correspond to features, and columns to examples.
# ind.pos : integer vector with the indices of the positive examples. Indices correspond to the columns of C.
# ind.markers : integer vector with the indices of the markers to be assessed. Indices correspond to the rows of W.
# score : Local score function to be considered. It may be one of the following:
#                 - eav.score (default)
#                 - NN.score
#                 - KNN.score 
#                 - WSLD score
# 		          - tot.score
# 		          - diff.score
#		          - dnorm.score
# kernel : kernel method (def. rw.kernel). 
# a : kernel parameter (def. 2).
# p : number of steps of the RW kernel (def. 1). It is meaningful only with rw.kernel. If 0, it is ignored. 
# probs: probabilities relative to the quantiles to be computed (def. 0.1 to 1 by 0.1 steps). Minimum values of W (usually zeros) are not considered.
#        Note: the numeric vector must be monotonically crescent. If intloo=FALSE, probs represents the quantile that will be used to filter the network.
# rep : number of the repetitions of the resamplings (def: 1000).
# intloo : logical. If TRUE (default) an internal loo is performed to find the optimal threshold, otherwise no internal loo is performed and the probs value is used to compute the quantile for filtering. 
# seed : init seed for the random resamplings (def. NULL), If NULL no initialization is performed.
# ... : optional arguments for the function score. 
# Output:
# a list with the following elements:
#   - marker.p.value : the computed p-value for the assessed marker set.
#   - auc : value of the computed AUC by loo.
#   - scores : vector of the scores computed for each example.
markerset.assessment.loo <-function(C, ind.pos, ind.markers, score=eav.score, kernel=rw.kernel, a=2, p=1, probs=seq(0.1, 1, 0.1), rep=1000, intloo=TRUE, seed=NULL, ...) {

 set.seed(seed);
 m <- ncol(C); # number of examples
 n <- nrow(C); # number of features
 
 # pnet computed on the true labels
 W <- C[ind.markers,];
 K <- cor(W);
 K[K<0 | is.na(K)] <- 0;
 K <- kernel(K, a);
 if (p>1)
    K <- p.step.rw.kernel(K, p=p); 
 res <- pnet.loo(K, ind.pos, score=score, probs=probs,  intloo=intloo, ...);
 s <- res$s;
 names(s) <- rownames(K); 
 q <- mean(res$q);
 target <- numeric(m);
 target[ind.pos] <- 1;
 auc <- AUC.single(s, target);
 
 # random resampling of markers
 marker.p.value <- 0.0;
 n.markers <- length(ind.markers);
 #x <- setdiff(1:n,ind.markers);
 x <- 1:n
 for (i in 1:rep) {
   ind.markers.random <- sample(x, n.markers);
   W <- C[ind.markers.random,];
   K <- cor(W);
   K[K<0 | is.na(K)] <- 0;
   K <- kernel(K, a);
   if (p>1)
     K <- p.step.rw.kernel(K, p=p);   
   res <- pnet.loo(K, ind.pos, score=score, probs=q,  intloo=FALSE, ...);
   auc.random <- AUC.single(res$s, target);
   if (auc.random >= auc)
     marker.p.value <- marker.p.value + 1;
 }
 return(list(marker.p.value=marker.p.value/rep, auc=auc, scores=s));
}

# Function that implements a non parametric test to assess the statistical significance of a given biomarker set through cross validation.
# The test is based on multiple random resamplings of the biomarker set, and a p.value is estimated by comparing the AUC 
# achieved with the assessed biomarker set with those obtained by random resampling.
# The threshold to filter the net for the external cross validation is obtained by averaging across the thresholds estimated at each internal loo.
# Note that if intloo=FALSE no internal loo is executed and the net is filtered according to the quantile of the argument probs.
# Input:
# C : a matrix representing examples and the associated features. Rows correspond to features, and columns to examples.
# ind.pos : integer vector with the indices of the positive examples. Indices correspond to the columns of C.
# ind.markers : integer vector with the indices of the markers to be assessed. Indices correspond to the rows of W.
# kk : number of folds (def. 5).
# score : Local score function to be considered. It may be one of the following:
#                 - eav.score (default)
#                 - NN.score
#                 - KNN.score 
#                 - WSLD score
# 		          - tot.score
# 		          - diff.score
#		          - dnorm.score
# kernel : kernel method (def. rw.kernel).
# a : kernel parameter (def. 2).
# p : number of steps of the RW kernel (def. 1). It is meaningful only with rw.kernel. If 0, it is ignored. 
# probs: probabilities relative to the quantiles to be computed (def. 0.1 to 1 by 0.1 steps). Minimum values of W (usually zeros) are not considered.
#        Note: the numeric vector must be monotonically crescent. If intloo=FALSE, probs represents the quantile that will be used to filter the network.
# rep : number of the repetitions of the resamplings (def: 1000).
# intloo : logical. If TRUE (default) an internal loo is performed to find the optimal threshold, otherwise no internal loo is performed and the probs value is used to compute the quantile for filtering. 
# seed : init seed for the random resamplings (def. NULL), If NULL no initialization is performed.
# ... : optional arguments for the function score. 
# Output:
# a list with the following elements:
#   - marker.p.value : the computed p-value for the assessed marker set.
#   - auc : value of the computed AUC by cross validation.
#   - scores : vector of the scores computed for each example.
markerset.assessment.cv <-function(C, ind.pos, ind.markers, kk=5, score=eav.score, kernel=rw.kernel, a=2, p=1, probs=seq(0.1, 1, 0.1), rep=1000, intloo=TRUE, seed=NULL, ...) {

 set.seed(seed);
 m <- ncol(C); # number of examples
 n <- nrow(C); # number of features
 
 # pnet computed on the true labels
 W <- C[ind.markers,];
 K <- cor(W);
 K[K<0] <- 0;
 K <- kernel(K, a);
 if (p>1)
    K <- p.step.rw.kernel(K, p=p); 
 res <- pnet.cv(K, ind.pos, kk=kk, score=score, probs=probs,  intloo=intloo, seed=seed, ...);
 s <- res$s;
 names(s) <- rownames(K); 
 q <- mean(res$q);
 target <- numeric(m);
 target[ind.pos] <- 1;
 auc <- AUC.single(s, target);
 
 # random resampling of markers
 marker.p.value <- 0.0;
 n.markers <- length(ind.markers);
 #x <- setdiff(1:n,ind.markers);
 x <- 1:n
 for (i in 1:rep) {
   ind.markers.random <- sample(x, n.markers);
   W <- C[ind.markers.random,];
   K <- cor(W);
   K[K<0] <- 0;
   K <- kernel(K, a);
   if (p>1)
     K <- p.step.rw.kernel(K, p=p);   
   res <- pnet.cv(K, ind.pos, kk=kk, score=score, probs=q,  intloo=FALSE, ...);
   auc.random <- AUC.single(res$s, target);
   if (auc.random >= auc)
     marker.p.value <- marker.p.value + 1;
 }
 return(list(marker.p.value=marker.p.value/rep, auc=auc, scores=s));
}


# Function that implements a non parametric test to assess the statistical significance of a given biomarker set through held out.
# The test is based on multiple random resamplings of the biomarker set, and a p.value is estimated by comparing the AUC 
# achieved with the assessed biomarker set with those obtained by random resampling.
# The threshold to filter the net for the test examples is obtained by loo on the training set.
# Note that if intloo=FALSE no internal loo is executed and the net is filtered according to the quantile of the argument probs.
# Input:
# C : a matrix representing examples and the associated features. Rows correspond to features, and columns to examples.
# ind.pos : integer vector with the indices of the positive examples. Indices correspond to the columns of C.
# ind.markers : integer vector with the indices of the markers to be assessed. Indices correspond to the rows of W.
# test : ratio of the test data (def. 0.3). I.e. if test is 0.3 of the available data, the train set will be the remaining 0.7.
#         Alternatively may be a vector with the indices of the test data. Indices corrrespond to the rows of W.
# score : Local score function to be considered. It may be one of the following:
#                 - eav.score (default)
#                 - NN.score
#                 - KNN.score 
#                 - WSLD score
# 		          - tot.score
# 		          - diff.score
#		          - dnorm.score
# kernel : kernel method (def. rw.kernel).
# a : kernel parameter (def. 2).
# p : number of steps of the RW kernel (def. 1). It is meaningful only with rw.kernel. If 0, it is ignored.
# probs: probabilities relative to the quantiles to be computed (def. 0.1 to 1 by 0.1 steps). Minimum values of W (usually zeros) are not considered.
#        Note: the numeric vector must be monotonically crescent. If intloo=FALSE, probs represents the quantile that will be used to filter the network.
# rep : number of the repetitions of the resamplings (def: 1000).
# intloo : logical. If TRUE (default) an internal loo is performed to find the optimal threshold, otherwise no internal loo is performed and the probs value is used to compute the quantile for filtering.
# seed : init seed for the random resamplings (def. NULL), If NULL no initialization is performed.
# ... : optional arguments for the function score.
# Output:
# a list with the following elements:
#   - marker.p.value : the computed p-value for the assessed marker set.
#   - auc : value of the computed AUC on the test data.
#   - scores : a vector of the scores (including both training and test data).
#   - test : indices of the elements of the test set (indices of the vector scores).
markerset.assessment.heldout <-function(C, ind.pos, ind.markers, test=0.3, score=eav.score, kernel=rw.kernel, a=2, p=1, probs=seq(0.1, 1, 0.1), rep=1000, intloo=TRUE, seed=NULL, ...) {

 set.seed(seed);
 m <- ncol(C); # number of examples
 n <- nrow(C); # number of features
 
 # pnet computed on the true labels
 W <- C[ind.markers,];
 K <- cor(W);
 K[K<0] <- 0;
 K <- kernel(K, a);
 if (p>1)
    K <- p.step.rw.kernel(K, p=p); 
	 
 # pnet computed on the true labels
 res <- pnet.heldout(K, ind.pos=ind.pos, test=test, score=score, probs=probs, intloo=intloo, seed=seed, ...);
 s <- res$s;
 names(s) <- rownames(K); 
 q <- res$q;
 target <- numeric(m);
 target[ind.pos] <- 1;
 ind.test <- res$test;
 s.test<-s[ind.test];
 target.ind.test<-target[ind.test];
 auc <- AUC.single(s.test, target.ind.test);
 
 # random resampling of markers
 marker.p.value <- 0.0;
 n.markers <- length(ind.markers);
 #x <- setdiff(1:n,ind.markers);
 x <- 1:n;
 for (i in 1:rep) {
   ind.markers.random <- sample(x, n.markers);
   W <- C[ind.markers.random,];
   K <- cor(W);
   K[K<0] <- 0;
   K <- kernel(K, a);
   if (p>1)
     K <- p.step.rw.kernel(K, p=p);   
   res <- pnet.heldout(K, ind.pos=ind.pos, test=ind.test, score=score, probs=q, intloo=FALSE, ...);
   s.test<-res$s[ind.test];
   auc.random <- AUC.single(s.test, target.ind.test);
   if (auc.random >= auc)
     marker.p.value <- marker.p.value + 1;
 }
 return(list(marker.p.value=marker.p.value/rep, auc=auc, scores=s, test=ind.test));
}


####################################################################
# Methods implementing:
# classifiers based on P-Net
####################################################################


# Method that implements the classifier based on the P-Net algorithm using held-out.
# A loo on the training set is performed to select the optimal threshold.
# The threshold to filter the net for the external held-out is obtained by averaging across the thresholds estimated at each internal loo
# for each fold.
# Note that if intloo=FALSE no internal loo is executed and the net is filtered according to the quantile of the argument probs.
# The classifier is constructed by finding the "optimal" score theshold: examples whose scores are above the selected threshold are positive otherwise negative.
# The score theshold is "optimal" w.r.t. a specific metric (F-score or accuracy). To assure an unbiased estimate of the score threshold, its assessment is performed on the train. Several quantiles are attempted and the best one used is used to classify the examples of the test set.
# Input: 
# W : a symmetric matrix representing similarity between examples. It can be a kernel or a correlation matrix, or any other matrix representing meaningful similarities between examples.
# ind.pos : integer vector with the indices of the positive examples. Indices correspond to the row (column) of W.
# test : ratio of the test data (def. 0.3). I.e. if test is 0.3 of the available data, the train set will be the remaining 0.7.
#         Alternatively may be a vector with the indices of the test data. Indices corrrespond to the rows of W.
# score : Local score function to be considered. It may be one of the following:
#                 - eav.score (default)
#                 - NN.score
#                 - KNN.score 
#                 - WSLD score
# 		          - tot.score
# 		          - diff.score
#		          - dnorm.score
# probs : probabilities relative to the quantiles to be computed (def. 0.01 to 1 by 0.01 steps). Minimum values of W (usually zeros) are not considered.
#        Note: the numeric vector must be monotonically crescent. If intloo=FALSE, probs represents the quantile that will be used to filter the network.
# score.probs : probabilities relative to the quantiles to be tested (def. 0.1 to 1 by 0.1 steps) to find the optimal score threshold for the classification (def. 0.01 to 1 by 0.01 steps).
# intloo : logical. If TRUE (default) an internal loo is performed to find the optimal threshold, otherwise no internal loo is performed and the probs value is used to compute the quantile for filtering. 
# opt.fun : function. Function implementing the metric to choice the optimal score threshold for classification. The F-score (compute.F) is the default.
#           Available functions:
#                 - compute.F (default) : F-score
#                 - compute.acc  : accuracy
#           N.B.: any function having two arguments representing the vector of predicted and true labels can be in principle used.
# seed : init seed for the random construction of the folds (def. 1). 
# ... : optional arguments for the function score. 
# Output:
# a list with the following elements:
#   - labels : a vector with the predicted labels (including both training and test data): 1 stands for positive and 0 for negative.
#   - s : a vector of the scores (including both training and test data).
#   - test : indices of the elements of the test set (indices of the vector s).
#   - q : the quantile  selected on the training set to filter the net.
#   - q.score : the quantile  score selected on the training set to classify the test set.
#   - F : the F-score computed on the test set.
#   - acc : the accuracy computed on the test set.
setGeneric("pnet.class.heldout", 
                 function(W, ind.pos, test=0.3, score=eav.score, probs=seq(0.1, 1, 0.1), score.probs=seq(0.01, 1, 0.01),
				          intloo=TRUE, opt.fun=compute.F, seed=NULL, ...) standardGeneric("pnet.class.heldout"));
				 
setMethod("pnet.class.heldout", signature(W="matrix"),
  function(W, ind.pos, test=0.3, score=eav.score, probs=seq(0.1, 1, 0.1), score.probs=seq(0.01, 1, 0.01), intloo=TRUE, opt.fun=compute.F, seed=NULL, ...) {
  
    res <- pnet.heldout(W, ind.pos, test=test, score=score, probs=probs, intloo=intloo, seed=seed, ...);
	m <- nrow(W);
	scores <- res$s;
	ind.test <- res$test;
	q <- res$q;

	labels <- integer(m);
	labels[ind.pos] <- 1;
	ind.train <- setdiff(1:m, ind.test);
	labels.train <- labels[ind.train];
	scores.train <- scores[ind.train];
	q.score=0;
	opt.thresh<-0;
	opt.F=0;
	for (qq in score.probs) {
	  thresh <- quantile(scores.train, probs=qq);
	  pred.train <- ifelse(scores.train>=thresh, 1, 0);
	  F <- opt.fun(pred.train,labels.train);    
	  if (F >= opt.F) {
	   opt.F <- F;
	   q.score <- qq;
	   opt.thresh <- thresh;
	  }	
	}
	# N.B.: the threshold is computed using the quantile computed on the train set but applied on the test set
	#       this should improve the sensitivity
	thresh <- quantile(scores[ind.test], probs=q.score);
	#thresh <- quantile(scores, probs=q.score);
	pred.labels <- ifelse(scores>=thresh, 1, 0);
	names(pred.labels)<-rownames(W);
	
	F<-compute.F(pred.labels[ind.test],labels[ind.test]);
	acc<-compute.acc(pred.labels[ind.test],labels[ind.test]);
	
    return (list(labels=pred.labels, s=scores, test=ind.test, q=q, q.score=q.score, F=F, acc=acc)); 
})


# Method that implements the classifier based on the P-Net algorithm using cross-validation.
# An external CV is performed to evaluate the performance of the method, and an internal loo is performed to select the optimal threshold.
# The threshold to filter the net for the external cv is obtained by averaging across the thresholds estimated at each internal loo
# for each fold.
# Note that if intloo=FALSE no internal loo is executed and the net is filtered according to the quantile of the argument probs.
# The classifier is constructed by finding the "optimal" score theshold: examples whose scores are above the selected threshold are positive otherwise negative.
# The score theshold is "optimal" w.r.t. a specific metric (F-score or accuracy). To assure an unbiased estimate of the score threshold, its assessment is performed on the train sets: several quantiles are tested and the best one is used to classify the examples of the corresponding test set.
# Input: 
# W : a symmetric matrix representing similarity between examples. It can be a kernel or a correlation matrix, or any other matrix representing meaningful similarities between examples. 
# ind.pos : integer vector with the indices of the positive examples. Indices correspond to the row (column) of W.
# kk : number of folds (def. 5).
# score : Local score function to be considered. It may be one of the following:
#                 - eav.score (default)
#                 - NN.score
#                 - KNN.score 
#                 - WSLD score
# 		          - tot.score
# 		          - diff.score
#		          - dnorm.score
# probs : probabilities relative to the quantiles to be computed (def. 0.1 to 1 by 0.1 steps). Minimum values of W (usually zeros) are not considered.
#        Note: the numeric vector must be monotonically crescent. If intloo=FALSE, probs represents the quantile that will be used to filter the network.
# score.probs: probabilities relative to the quantiles to be tested (def. 0.01 to 1 by 0.01 steps) to find the optimal score threshold for the classification (def. 0.01 to 1 by 0.01 steps).
# intloo : logical. If TRUE (default) an internal loo is performed to find the optimal threshold, otherwise no internal loo is performed and the probs value is used to compute the quantile for filtering. 
# opt.fun : function. Function implementing the metric to choice the optimal score threshold for classification. The F-score (compute.F) is the default.
#           Available functions:
#                 - compute.F (default) : F-score
#                 - compute.acc  : accuracy
#           N.B.: any function having two arguments representing the vector of predicted and true labels can be in principle used.
# seed : init seed for the random construction of the folds (def. 1).
# ... : optional arguments for the function score. 
# Output:
# a list with the following elements:
#   - labels : vector of labels computed by k-fold cross validation. 
#   - s : a vector with scores computed by k-fold cross validation.
#   - q : the vector of kk quantiles corresponding to the quantiles selected at each fold to filter the net.
#   - q.score : the vector of kk quantile scores selected on the training set to classify the test set for each fold. 
#   - F : the F-score computed by k-fold cross validation. 
#   - acc : the accuracy computed by k-fold cross validation. 

setGeneric("pnet.class.cv", 
                 function(W, ind.pos, kk=5, score=eav.score, probs=seq(0.1, 1, 0.1), score.probs=seq(0.01, 1, 0.01), intloo=TRUE, opt.fun=compute.F, seed=1, ...) standardGeneric("pnet.class.cv"));
				 
setMethod("pnet.class.cv", signature(W="matrix"), 
  function(W, ind.pos, kk=5, score=eav.score, probs=seq(0.1, 1, 0.1), score.probs=seq(0.01, 1, 0.01), intloo=TRUE, opt.fun=compute.F, seed=1, ...) {
    
	score.name <- return.name(score);
	# filtering of the net
    res <- pnet.cv(W, ind.pos, kk=5, score=score, probs=probs, intloo=intloo, seed=seed, ...);
	
	m <- nrow(W);
	scores <- res$s;
	q <- res$q;
	q.score <- numeric(kk);
	x <- 1:m;
	thresh <- quantile(W,probs=mean(q));  # reconstruction of the W matrix	
	W[W<thresh] <- 0;
	folds <- do.stratified.cv.data(x, ind.pos, k=kk, seed=seed); # reconstruction of the folds
	labels <- pred.labels <- integer(m);
	names(pred.labels) <- names(scores) <-rownames(W);
	labels[ind.pos] <- 1;
	
	# cycling on the kk folds: score thresh (quantiles) are computed on the training data and used to classify on the test set
	for (j in 1:kk) {
	  ind.test <- c(folds$fold.positives[[j]], folds$fold.non.positives[[j]]);
	  ind.train <- x[-ind.test];
	  ind.pos.train <- setdiff(ind.pos,ind.test);	  
	  labels.train <- labels[ind.train];	  
          if (score.name == "KNN.score")  				  
 	  	scores.train <- score(W, ind.train, ind.pos.train, k=k, norm=FALSE, ...)   
    	  else if (score.name == "WSLD.score")		  
 	  	scores.train <- score(W, ind.train, ind.pos.train, d=d, norm=FALSE, ...) 
    	  else if (score.name == "tot.score" || score.name == "diff.score" || score.name == "dnorm.score") 
          	scores.train <- score(W, ind.train, ind.pos.train) 
    	  else										  
 	  	scores.train <- score(W, ind.train, ind.pos.train, norm=FALSE, ...); 
	  q.score[j] <- 0;
	  opt.F <- 0;
	  for (qq in score.probs) {
	     thresh <- quantile(scores.train, probs=qq);
	     pred.train <- ifelse(scores.train>=thresh, 1, 0);
	     F <- opt.fun(pred.train,labels.train);    
	     if (F >= opt.F) {
	       opt.F <- F;
	       q.score[j] <- qq;
	     } 
	  }
	  # N.B.: the threshold is computed using the quantile computed on the train set but applied on the test set
	  # this should improve the sensitivity
	  thresh <- quantile(scores[ind.test], probs=q.score[j]);	
	  pred.labels[ind.test] <- ifelse(scores[ind.test]>=thresh, 1, 0);
	}
	
    F<-compute.F(pred.labels,labels);
	acc<-compute.acc(pred.labels,labels);
	
    return (list(labels=pred.labels, s=scores, q=q, q.score=q.score, F=F, acc=acc)); 
	
})


# Method that implements the classifier based on the P-Net algorithm using loo.
# An external loo is performed to evaluate the performance of the method, and an internal loo is performed to select the optimal threshold.
# The threshold to filter the net for the external loo is obtained by averaging across the thresholds estimated at each internal loo
# for each fold.
# Note that if intloo=FALSE no internal loo is executed and the net is filtered according to the quantile of the argument probs.
# The classifier is constructed by finding the "optimal" score theshold: examples whose scores are above the selected threshold are positive otherwise negative.
# The score threshold is "optimal" w.r.t. a specific metric (F-score or accuracy). To assure an unbiased estimate of the score threshold, its assessment is performed on the train sets (data without the left out example): several quantiles are tested and the best one is used to classify the examples of the corresponding left out example.
# Input: 
# W : a symmetric matrix representing similarity between examples. It can be a kernel or a correlation matrix, or any other matrix representing meaningful similarities between examples. 
# ind.pos : integer vector with the indices of the positive examples. Indices correspond to the row (column) of W.
# score : Local score function to be considered. It may be one of the following:
#                 - eav.score (default)
#                 - NN.score
#                 - KNN.score 
#                 - WSLD score
# 		          - tot.score
# 		          - diff.score
#		          - dnorm.score
# probs: probabilities relative to the quantiles to be computed (def. 0.1 to 1 by 0.1 steps). Minimum values of W (usually zeros) are not considered.
#        Note: the numeric vector must be monotonically crescent. If intloo=FALSE, probs represents the quantile that will be used to filter the network.
# score.probs: probabilities relative to the quantiles to be tested (def. 0.01 to 1 by 0.01 steps) to find the optimal score threshold for the classification (def. 0.01 to 1 by 0.01 steps).
# intloo : logical. If TRUE (default) an internal loo is performed to find the optimal threshold, otherwise no internal loo is performed and the probs value is used to compute the quantile for filtering. 
# opt.fun : function. Function implementing the metric to choice the optimal score threshold for classification. The F-score (compute.F) is the default.
#           Available functions:
#                 - compute.F (default) : F-score
#                 - compute.acc  : accuracy
#           N.B.: any function having two arguments representing the vector of predicted and true labels can be in principle used.
# ... : optional arguments for the function score. 
# Output:
# a list with the following elements:
#   - labels : vector of labels computed by loo.
#   - s : a vector with scores computed by loo.
#   - q : the vector of quantiles selected for each example to filter the net.
#   - q.score : the vector of quantile scores selected on the training set to classify the left out example.
#   - F : the F-score computed by loo.
#   - acc : the accuracy computed by loo. 

setGeneric("pnet.class.loo", 
                 function(W, ind.pos, score=eav.score, probs=seq(0.1, 1, 0.1), score.probs=seq(0.01, 1, 0.01), intloo=TRUE, opt.fun=compute.F, ...) standardGeneric("pnet.class.loo"));
				 
setMethod("pnet.class.loo", signature(W="matrix"), 
  function(W, ind.pos, score=eav.score, probs=seq(0.1, 1, 0.1), score.probs=seq(0.01, 1, 0.01), intloo=TRUE, opt.fun=compute.F, ...) {
    
	score.name <- return.name(score);
	# filtering of the net
    res <- pnet.loo(W, ind.pos,  score=score, probs=probs, intloo=intloo,  ...);
	
	m <- nrow(W);
	scores <- res$s;
	q <- mean(res$q);
	q.score <- numeric(m);
	x <- 1:m;
	labels <- pred.labels <- integer(m);
	names(pred.labels) <- names(scores) <-rownames(W);
	labels[ind.pos] <- 1;
    thresh <- quantile(W,probs=q);  # reconstruction of the W matrix	
	W[W<thresh] <- 0;
	
	# cycling on the kk folds: score thresh (quantiles) are computed on the training data and used to classify on the test set
	for (j in 1:m) {
          if (score.name == "KNN.score")  				  
 	  	scores.train <- score(W, x[-j], setdiff(ind.pos,j), k=k, norm=FALSE, ...)   
    	  else if (score.name == "WSLD.score")		  
 	  	scores.train <- score(W, x[-j], setdiff(ind.pos,j), d=d, norm=FALSE, ...) 
    	  else if (score.name == "tot.score" || score.name == "diff.score" || score.name == "dnorm.score") 
          	scores.train <- score(W, x[-j], setdiff(ind.pos,j)) 
    	  else										  
 	  	scores.train <- score(W, x[-j], setdiff(ind.pos,j), norm=FALSE, ...); 
	  labels.train <- labels[-j];	  
	  q.score[j] <- 0;
	  opt.F <- 0;
	  for (qq in score.probs) {
	     thresh <- quantile(scores.train, probs=qq);
	     pred.train <- ifelse(scores.train>=thresh, 1, 0);
	     F <- opt.fun(pred.train,labels.train);    
	     if (F >= opt.F) {
	       opt.F <- F;
	       q.score[j] <- qq;
	     } 
	  }
	  thresh <- quantile(scores, probs=q.score[j]);	
	  pred.labels[j] <- ifelse(scores[j]>=thresh, 1, 0);
	}
	
    F<-compute.F(pred.labels,labels);
	acc<-compute.acc(pred.labels,labels);
	
    return (list(labels=pred.labels, s=scores, q=q, q.score=q.score, F=F, acc=acc)); 
	
})

###########################################################################
# Functions implementing the visual representation of the results as graphs
###########################################################################


# Function to plot results obtained through the P-Net algorithm.
# This function plots a graph, representing patients as nodes and similarities between patients through edges.
# Positive patients, i.e. patients having the phenotype of interest are represented as squares, while the other patients are represented with circles.
# The scores obtained by cross validation by the algoritm are represented through shaded colours, using a base color (default:red).
# Patients with the highest scores are represented through more intense colours, while patients having the lowest scores 
# with light shades of color.
# Input:
#   A : a square symmetric matrix representing an undirected graph. It is an adjiacency matrix whose values can be 1 (edge) or 0 (no edge).
#   ind.pos: a vector with the indices of a priori known positive patients.
#   scores : vector of the scores computed through P-Net. 
#   colour: colour of the nodes to be highlighted (def: red). Allowed colors are: red, green, blue, gray.
#   file.name : name of the postscript file of the constructed graph. If "" (def.) the graph is printed into a window.
#   mode : mode of layout. It can be "dot" (default), "twopi", "neato".
#   fontsize: a real value. Size of the font used for the label of the nodes (def:12).
#   magn : a real value. It is the "magnification" of the  nodes (def.1, that is no magnification).
# Output:
#   The plot of the graph is printed into a window or as a postscript file.
plot.pnet.graph <- function (A, ind.pos, scores, colour="red", file.name="", mode="dot", fontsize=12, magn=1) {
  library(graph);
  library(Rgraphviz);
  library("RColorBrewer");
  
  if (ncol(A) != nrow(A))
    stop("Adjagency matrix must be square");
  
  n <- nrow(A); 
  rownames(A)<-colnames(A)<-names(scores)<-1:n;  # this is to set the node labels to integers
  # Making graph
  gw <- new("graphAM", adjMat=A,values=list(weight=1));
  
  # setting global rendering attributes of the graph
  defAttrs <- getDefaultAttrs();
  nheight <- as.numeric(defAttrs$node$width);
  nwidth <- as.numeric(defAttrs$node$width);
  attrs=list(node=list(fillcolor="white", fontsize=fontsize, height=nheight*magn, width=nwidth*magn), 
             edge=list(color="gray"), graph=list(rankdir="TB"));
  
  # setting per node attributes
  nAttrs <- list();
  if (colour == "red")
     colors <- colorRampPalette(brewer.pal(9, "Reds"))(255)   
  else if (colour == "blue")
     colors <- colorRampPalette(brewer.pal(9, "Blues"))(255)
  else if (colour == "green")
     colors <- colorRampPalette(brewer.pal(9, "Greens"))(255)
  else if (colour == "gray")
     colors <- colorRampPalette(brewer.pal(9, "Greys"))(255)
  else if (colour == "purple")
     colors <- colorRampPalette(brewer.pal(9, "Purples"))(255)
  else
     stop("plot.pnet.graph; color not supported");
  ind.colors<-integer(n);
  max.score <- max(scores);
  for (i in 1:n) 
    ind.colors[i] <- as.integer(1 + round((scores[i]*(255-1))/max.score)); 
  node.fillcolor <- colors[ind.colors];
  names(node.fillcolor) <- nodes(gw);
  nAttrs$fillcolor <- node.fillcolor;
  
  nodes.shape <- rep("circle",n);
  nodes.shape[ind.pos] <- "box";
  names(nodes.shape) = nodes(gw);
  nAttrs$shape <- nodes.shape;
  
  # plot
  if (file.name=="") {
    x11(); 
    plot(gw, nodeAttrs = nAttrs, attrs=attrs, y=mode);
  } else {
    postscript(file.name, width=12, height=8, horizontal=TRUE);
    plot(gw, nodeAttrs = nAttrs, attrs=attrs, y=mode);
    dev.off();
  }  
}


##########################################################
#Useful functions not present in RANKS 1.0
##########################################################

# Function that return the string name of the function included in a variable.
# Input:
# fun : variable storing a function.
# Output:
# a string with the name of the function.
return.name <- function(fun) {
   x<-showMethods(fun, printTo=FALSE);
   fun.name <- strsplit(x[1],split=" ")[[1]][2];
   return(fun.name)
 }

# Function to compute the best score threshold.
# Input:
# s : scores.
# labels : true labels (1 positive, 0 negative).
# opt.fun : function. Function implementing the metric to choose the optimal score threshold for classification. The accuracy (compute.acc) is the default.
#           Available functions:
#                 - compute.F  : F-score
#                 - compute.acc  : accuracy
#           N.B.: any function having two arguments representing the vector of predicted and true labels can be in principle used.
# score.probs: probabilities relative to the quantiles to be tested (def. 0.1 to 1 by 0.1 steps) to find the optimal score threshold for the classification (def. 0.01 to 1 by 0.01 steps).
# Output:
# the best quantile for score thresholding 
find.best.quantile.score <- function (s, labels, opt.fun=compute.acc, score.probs=seq(0.01, 1, 0.01)) {
   opt.q.score=0;
   opt.F=0;
   for (qq in score.probs) {
     thresh <- quantile(s, probs=qq);
     pred <- ifelse(s>=thresh, 1, 0);
     F <- opt.fun(pred,labels);    
     if (F >= opt.F) {
      opt.F <- F;
      opt.q.score <- qq;
     } 
   }
   return(opt.q.score);
}

# Function to generate data for the stratified hold.out. 
# Input:
# examples : indices of the examples (a vector of integer).
# positives: vector of integer. Indices of the positive examples. The indices refer to the indices (values) of examples.
# train.ratio : ratio of the training data w.r.t. the overall number of available data (def = 0.7).
# seed : seed of the random generator (def=0). If is set to 0 no initialitation is performed.
# Ouptut:
# a list with 4 components
#   - train.pos : vector with the indices of the positive elements of the training set.
#   - train.neg : vector with the indices of the negative elements of the training set.
#   - test.pos : vector with the indices of the positive elements of the test set.
#   - test.neg : vector with the indices of the negative elements of the test set.
# N.B.: in all the components of the list, the elements are those included in examples (i.e. values of the vector examples) 
do.stratified.hold.out.data <- function(examples, positives, train.ratio=0.7, seed=0) {
  if (seed!=0)
    set.seed(seed);
  res <- list();
  for (i in 1:4)
    res[[i]] <- integer(0);
  names(res) <- c("train.pos", "train.neg", "test.pos", "test.neg");
  
  positives <- sample(positives);
  n.positives <- length(positives);
  negatives <- setdiff(examples,positives);
  negatives <- sample(negatives);
  n.negatives <- length(negatives);
  if ((n.positives==0) || (n.negatives==0))
    stop("do.stratified.hold.out.data: no positives or negatives in examples!");
  
  n.positives.train <- round(n.positives * train.ratio);
  n.positives.test <- n.positives - n.positives.train;
  n.negatives.train <- round(n.negatives * train.ratio);
  n.negatives.test <- n.negatives - n.negatives.train;
  
  res$train.pos <- sort(positives[1:n.positives.train]);
  if (n.positives.train == n.positives)
     res$test.pos <- integer(0)
  else 
     res$test.pos <- sort(positives[(n.positives.train+1):n.positives]);

  res$train.neg <- sort(negatives[1:n.negatives.train]);
  if (n.negatives.train == n.negatives)
     res$test.neg <- integer(0)
  else 
     res$test.neg <- sort(negatives[(n.negatives.train+1):n.negatives]); 
  return(res);
}

# Function to compute the F-measure, precision and recall for a single class.
# Input:
# pred : factor of the predicted labels.
# labels : factor of the true labels.
# Note that 0 level stands for negative and 1 for positive.
# In general the first level is negative and the second positive.
# Output:
# A vector with F, precision and recall.
compute.F.prec.rec <- function(pred,labels) {      
     if (length(pred)!=length(labels))
         stop("compute.F: lengths of true and predicted labels do not match.");
     neg.ex <- which(labels <= 0);	
	 pos.ex <- which(labels > 0);
	 TP <- sum(pred[pos.ex] > 0);
	 FN <- sum(pred[pos.ex] <= 0);	
	 TN <- sum(pred[neg.ex] <= 0);
	 FP <- sum(pred[neg.ex] > 0);	           
     if ((TP+FP) == 0)
       precision <- 0
     else 
       precision <- TP/(TP+FP);
     if ((TP+FN) == 0)
       recall <- 0
     else
       recall <- TP/(TP+FN);
     if ((precision+recall) == 0)
     	F <- 0
     else
     	F = 2 *(precision*recall) / (precision+recall);    
	 res <- c(F,precision,recall); 
	 names(res) <- c("F", "P", "R"); 
     return (res);
}

###########################################################
#Score functions not present in RANKS_1.0
###########################################################

#TOTAL SCORE
setGeneric("tot.score", 
                 function(RW, x, x.pos) standardGeneric("tot.score"));

# Method to compute the Total score for a set of vertices.
# score(x) = \frac{ \sum_{x_i \in x.pos} K(x,x_i)}{\sum_{x_i \in x.pos} K(x,x_i) + \sum_{x_i \in x[-x.pos]} K(x,x_i)}}
# Input:
# RW : matrix. It must be a kernel matrix or a symmetric matrix representing the similarity
#              between pairs of nodes.
# x : vector of integer. Indices corresponding to the elements of the RW matrix for which the score must be computed.
# x.pos : vector of integer. Indices of the positive elements of the RW matrix. 
# Output:
# vector of the total scores of the elements x. The names of the vector correspond to the indices of x.
setMethod("tot.score", signature(RW="matrix"),
  function(RW, x, x.pos) {
	if (nrow(RW) != ncol(RW))
            stop("second arg must be a square matrix");
            
        score <- (apply(as.matrix(RW[x,x.pos]),1,sum))/((apply(as.matrix(RW[x,x.pos]),1,sum))+(apply(as.matrix(RW[x,x[-x.pos]]),1,sum)));
	                                              
        names(score)<-x;  
        score[is.nan(score)] <- 0;
	return(score);
})


# Method to compute the Total score for a set of vertices. 
# Input:
# x : vector of integer. Indices corresponding to the elements of the RW matrix for which the score must be computed.
# RW : an object of the virtual class graph (hence including objects of class graphAM  and graphNEL from the package graph).
# x.pos : vector of integer. Indices of the positive elements of the RW matrix.
# Output:
# vector of the total scores of the elements x. The names of the vector correspond to the indices of x.
setMethod("tot.score", signature(RW="graph"),
  function(RW, x, x.pos) {
     RW <- as(RW, "matrix");
     return(tot.score(RW, x, x.pos));
})

#DIFFERENTIAL SCORE
setGeneric("diff.score", 
                 function(RW, x, x.pos) standardGeneric("diff.score"));

# Method to compute the Differential score for a set of vertices.
# score(x) = \sum_{x_i \in x.pos} K{x,x_i} - \sum_{x_i \in x[-x.pos]} K{x,x_i} 
# Input:
# RW : matrix. It must be a kernel matrix or a symmetric matrix representing the similarity
#              between pairs of nodes.
# x : vector of integer. Indices corresponding to the elements of the RW matrix for which the score must be computed. 
# x.pos : vector of integer. Indices of the positive elements of the RW matrix. 
# Output: 
# vector of the differential scores of the elements x. The names of the vector correspond to the indices of x. 
setMethod("diff.score", signature(RW="matrix"),
  function(RW, x, x.pos) {
	if (nrow(RW) != ncol(RW))
            stop("second arg must be a square matrix");
	
        score <- (apply(as.matrix(RW[x,x.pos]),1,sum)) - (apply(as.matrix(RW[x,x[-x.pos]]),1,sum));
	names(score)<-x;  
        score[is.nan(score)] <- 0;
	return(score);
})


# Method to compute the Differential score for a set of vertices. 
# Input:
# x : vector of integer. Indices corresponding to the elements of the RW matrix for which the score must be computed. 
# RW : an object of the virtual class graph (hence including objects of class graphAM and graphNEL from the package graph).
# x.pos : vector of integer. Indices of the positive elements of the RW matrix. 
# Output:
# vector of the differential scores of the elements x. The names of the vector correspond to the indices of x. 
setMethod("diff.score", signature(RW="graph"),
  function(RW, x, x.pos) {
     RW <- as(RW, "matrix");
     return(diff.score(RW, x, x.pos));
})

#DIFFERENTIAL NORMALIZED SCORE
setGeneric("dnorm.score", 
                 function(RW, x, x.pos) standardGeneric("dnorm.score"));

# Method to compute the Differential Normalized score for a set of vertices.
# score(x) = \frac{\sum_{x_i \in x.pos} K{x,x_i} - \sum_{x_i \in x[-x.pos]} K{x,x_i}} {\sum_{x_i \in x.pos} K{x,x_i} + \sum_{x_i \in x[-ind.pos]} K{x,x_i}}
# Input:
# RW : matrix. It must be a kernel matrix or a symmetric matrix representing the similarity
#              between pairs of nodes.
# x : vector of integer. Indices corresponding to the elements of the RW matrix for which the score must be computed.
# x.pos : vector of integer. Indices of the positive elements of the RW matrix. 
# Output:
# vector of the differential normalized scores of the elements x. The names of the vector correspond to the indices of x. 
setMethod("dnorm.score", signature(RW="matrix"),
  function(RW, x, x.pos) {
	if (nrow(RW) != ncol(RW))
            stop("second arg must be a square matrix");
        score <- ((apply(as.matrix(RW[x,x.pos]),1,sum))-(apply(as.matrix(RW[x,x[-x.pos]]),1,sum)))/((apply(as.matrix(RW[x,x.pos]),1,sum))+(apply(as.matrix(RW[x,x[-x.pos]]),1,sum)));
	names(score)<-x;  
        score[is.nan(score)] <- 0;
	return(score);
})


# Method to compute the Differential Normalized score for a set of vertices.
# Input:
# x : vector of integer. Indices corresponding to the elements of the RW matrix for which the score must be computed.
# RW : an object of the virtual class graph (hence including objects of class graphAM and graphNEL from the package graph).
# x.pos : vector of integer. Indices of the positive elements of the RW matrix.
# Output:
# vector of the differential normalized scores of the elements x. The names of the vector correspond to the indices of x. 
setMethod("dnorm.score", signature(RW="graph"),
  function(RW, x, x.pos) {
     RW <- as(RW, "matrix");
     return(dnorm.score(RW, x, x.pos));
})

