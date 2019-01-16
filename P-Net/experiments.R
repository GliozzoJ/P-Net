# High level functions and methods to perform the experiments/workflows 
# employed in the paper "Network modeling of patients’ biomolecular profiles 
# for clinical phenotype/outcome prediction" submitted to Scientific Reports - Nature.


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



###############################################################################
# Methods implementing the experimental workflows
###############################################################################

#Function "pnet.ttest.mccv": feature selection, through T-test, and examples ranking through multiple hold-out (Monte Carlo cross-validation).
# This function both selects the features using T-test and ranks examples using P-net and the selected features. The ranking is performed with respect to the test set, while feature selection and net-threshold selection are performed through efficient internal loo on the training set. When the set of features have been selected as well as the threshold to cut the edges, the optimal threshold for the scores is estimated by using the scores previously estimated on the training set.
# The hold-out can be repeated and the resulting AUC, accuracy, F-score, precision and recall on the test set are stored. 
# For each repetition a set of features (always of length = n.feat) is selected and stored.
# Input:
# M : matrix of examples data. Rows correspond to features (e.g. expression for a given gene), columns to examples (patients).
# ind.pos : integer vector with the indices of the positive examples. Indices correspond to the columns of M.
# train.ratio : ratio of the training set data w.r.t. the total number of examples (def. 0.7).
# rep : number of repetitions of the hold-out procedure (def. 1).
# sim : similarity function to be used between patients (def: cor). 
# cor.method : method for the choice of the correlation type (def. "pearson"). Other possibilities are "kendall" or "spearman".
# score: function. It must be a kernel-based score method:
#                 - KNN.score 
#                 - NN.score 
#                 - eav.score  (default)
#                 - WSLD.score
# 		          - tot.score
# 		          - diff.score
#		          - dnorm.score
# k : number of neighbours for the KNN.score (def.5). Meaningful only if the KNN.score is selected.
# d : parameter for the WSLD.score (def. 2). Meaningful only if the WSLD.score is selected.
# ker : list containing the kernel methods. Possible kernels:
#       - rw.kernel (default)
#       - identity.kernel
#       - linear.kernel
#       - gaussian.kernel
#       - laplacian.kernel
#       - cauchy.kernel
#       - inv.multiquadric.kernel
#       - poly.kernel
# p : number of steps of the rw kernel. It is meaningful and is used only with rw.kernel. If 0 (def) it is ignored.
# probs: probabilities relative to the quantiles to be computed (def. 0.1 to 1 by 0.1 steps) to optimally cut the edges. Minimum values of W (usually zeros) are not considered.
#        Note: the numeric vector must be monotonically crescent.
# opt.fun : function. Function implementing the metric to choice the optimal score threshold for classification. The accuracy (compute.acc) is the default.
#           Available functions:
#                 - compute.F  : F-score
#                 - compute.acc : accuracy (def.)
#           N.B.: any function having two arguments representing the vector of predicted and true labels can be in principle used.
# score.probs: probabilities relative to the quantiles to be tested (def. 0.1 to 1 by 0.1 steps) to find the optimal score threshold for the classification (def. 0.01 to 1 by 0.01 steps).
# seed : initialization of the random generator to take into account the randomness due to ties. If seed=0 no initialization (def: 1).
# n.feat: number of selected features based on T-test computed on the training set. 
# ... : parameters of the kernel ker and of the function score.
# Output: 
# a list with the following elements:
#   - avg.res : a vector with the following components: average AUC, accuracy, F-score, precision, recall across repetitions.
#   - st.dev : a vector with the standard deviation of the above measures.
#   - res.full : a matrix with a number of rows equal to the number of repetitions rep, and having as columns AUC, accuracy, F-score, precision, recall obtained at each repetition.
#   - sel.feat: a list where each component is an integer vector with the indices of the selected features. Indices refer to the rows of M. The length of the list is equal to rep.
#   - sel.quantile: a vector whose components are the quantile corresponding to the computed network edge threshold at each repetition.
#   - sel.score.quantile : a vector whose components are the quantile corresponding to the best score threshold at each repetition.

pnet.ttest.mccv <- function(M, ind.pos, train.ratio=0.7, rep=1, sim=cor, cor.method="pearson", score=eav.score, k=5, d=2, ker=list(rw.kernel=rw.kernel), p=0, probs=seq(0.1,1,0.1), opt.fun=compute.acc, score.probs=seq(0.01, 1, 0.01), seed=1, n.feat=1000, ...) {

 if (seed!=0)
   set.seed(seed);
 n <- nrow(M); # number of features
 n.ex <- ncol(M); #number of examples 
 
 avg.res <- st.dev <- numeric(5);
 names(avg.res) <- names(st.dev) <- c("AUC", "Acc", "F", "Prec", "Rec");
 res.full <- matrix(numeric(length(avg.res)*rep), nrow=rep);
 colnames(res.full) <- names(avg.res);
 sel.quantile <- sel.score.quantile <- numeric(rep);
 s.test <- vector ("list", rep);
 s.train <- vector("list", rep);

 sel.feat <- vector("list", rep); 
 x <- 1:n.ex;
 ker.name <- names(ker);
 score.name <- return.name(score);
 ind.pos <- sort(ind.pos);
 labels <- integer(n.ex);
 labels[ind.pos] <- 1;

 # repeated multiple hold-out
 for (i in 1:rep) {  
    f <- do.stratified.hold.out.data(x, ind.pos, train.ratio=train.ratio, seed=seed+i);
  
	# feature selection with T-test on the training data:
	ind.train <- sort(x[c(f$train.pos, f$train.neg)]);
	ind.pos.train <-f$train.pos;
        labels.train <- labels[ind.train];
        p.value <- integer(n);

        for (k in 1:n) #t-test applied only to the training set
             p.value[k] <- t.test( M[k,ind.train]~labels.train, mu=0, alternative="two.sided", conf.level=0.95, var.eq=F, paired=F)$p.value;

        sorted <- sort.int(p.value, index.return=TRUE); 
	sorted.index <- integer(n);
	sf <- sel.feat[[i]] <- sorted.index <- sorted[[2]][1:n.feat];
        
	# computing the kernel on the selected features:
	if (ker.name == "rw.kernel"  || ker.name == "identity.kernel") {    
      C <- sim(M[sf,], method=cor.method);  					    
	  C[C<0 | is.na(C)] <- 0;								    
	  K <- ker[[1]](C, ...); 											    
      if (p>1)  													    
    	 K <- p.step.rw.kernel(K,p); 								    
     } else  														    
        K <- ker[[1]](t(M[sf,]), ...);
 									    
	# optimization of the edge threshold by internal loo on the training set
	if (score.name == "KNN.score")
	    res <- optimize.thresh.by.loo(K, ind.pos=ind.pos.train, ind.test=ind.train, score=score, probs=probs, norm=FALSE, k=k, ...)
	else if (score.name == "WSLD.score")								   
	    res <- optimize.thresh.by.loo(K, ind.pos=ind.pos.train, ind.test=ind.train, score=score, probs=probs, norm=FALSE, d=d, ...)
	else if (score.name == "tot.score" || score.name == "diff.score" || score.name == "dnorm.score")			   
            res <- optimize.thresh.by.loo(K, ind.pos=ind.pos.train, ind.test=ind.train, score=score, probs=probs)
	else
            res <- optimize.thresh.by.loo(K, ind.pos=ind.pos.train, ind.test=ind.train, score=score, probs=probs, norm=FALSE, ...);
    sel.quantile[i] <- res$quantile;
																		   
    # net filtering for test set 																						   
    thresh <- quantile(K, probs=res$quantile);  		  													   
    K[K<thresh]<-0; 
	# Estimation of the scores on the training set and of the optimal score theshold by internal loo on the training set	
	K1 <- K;
	diag(K1) <- 0;  
	if (score.name == "KNN.score")  				  
 	    scores <- score (K1, x=x[ind.train], x.pos=ind.pos.train, k=k, norm=FALSE, ...)   
        else if (score.name == "WSLD.score")		  
 	    scores <- score (K1, x=x[ind.train], x.pos=ind.pos.train, d=d, norm=FALSE, ...) 
        else if (score.name == "tot.score" || score.name == "diff.score" || score.name == "dnorm.score") 
            scores <- score (K1, x=x[ind.train], x.pos=ind.pos.train) 
        else										  
 	    scores <- score (K1, x=x[ind.train], x.pos=ind.pos.train, norm=FALSE, ...);   
	labels.train <- labels[ind.train];
	best.quantile.score <- sel.score.quantile[i] <- find.best.quantile.score(scores, labels.train, opt.fun=opt.fun, score.probs=score.probs);

        s.train[[i]] <- scores;
	# Estimation of the scores on the test set 
	if (score.name == "KNN.score")  				  
 	    scores <- score (K, x=x[-ind.train], x.pos=ind.pos.train, k=k, norm=FALSE, ...)   
        else if (score.name == "WSLD.score")		  
 	    scores <- score (K, x=x[-ind.train], x.pos=ind.pos.train, d=d, norm=FALSE, ...) 
        else if (score.name == "tot.score" || score.name == "diff.score" || score.name == "dnorm.score") 
            scores <- score (K, x=x[-ind.train], x.pos=ind.pos.train)   
        else										 
 	    scores <- score (K, x=x[-ind.train], x.pos=ind.pos.train, norm=FALSE, ...); 

        s.test[[i]] <- scores;
	# Computing performances on the test set
	labels.test <- labels[-ind.train];
	res.full[i, "AUC"] <- AUC.single (scores, labels.test);
	thresh <- quantile(scores, probs=best.quantile.score);
	pred.labels <- ifelse(scores>=thresh, 1, 0);
	res.full[i, "Acc"]<-compute.acc(pred.labels, labels.test);
	res.full[i, 3:5] <- compute.F.prec.rec(pred.labels, labels.test);
	cat("End repetition ... ", i, "\n");
  }  # end for multiple repetitions
  avg.res <- apply(res.full,2,mean);
  st.dev <- apply(res.full,2,sd);
  return(list(avg.res=avg.res, st.dev=st.dev, res.full=res.full, sel.feat=sel.feat, sel.quantile=sel.quantile, sel.score.quantile=sel.score.quantile, s.test=s.test, s.train=s.train));				   

}


# Function "pnet.modtstat.repcv": Feature selection, through Moderated T-statistic Method, and examples ranking through multiple k-fold CV (cross-validation).
# This function both selects the features using moderated t-statistic method and ranks examples using P-Net and the selected features.
# The ranking is performed with respect to the test set, while feature selection and net-threshold selection are performed through efficient internal loo on the training set. 
# When the set of features have been selected as well as the threshold to cut the edges, the optimal threshold for the scores is estimated by using the scores previously estimated on the training set.
# The k-fold CV can be repeated and the resulting AUC, accuracy, F-score, precision and recall on the test set are stored. 
# For each repetition a set of features (always of length = n.feat) and the predicted labels are stored.
# Input:
# M : matrix of examples data. Rows correspond to features (e.g. expression for a given gene), columns to examples (patients).
# ind.pos : integer vector with the indices of the positive examples. Indices corresponding to the columns of M.
# kk : number of folds for the k-fold CV (def. 5). 
# rep : number of repetitions of the k-fold CV procedure (def. 1). 
# sim : similarity function to be used between patients (def: cor). 
# cor.method : method for the choice of the correlation type (def. "pearson"). Other possibilities are "kendall" or "spearman". 
# score: function. It must be a kernel-based score method:
#                 - KNN.score 
#                 - NN.score 
#                 - eav.score (default)
#                 - WSLD.score
# 		          - tot.score
# 		          - diff.score
#		          - dnorm.score
# ker : list containing the kernel method. Possible kernels:
#       - rw.kernel (default)
#       - identity.kernel
#       - linear.kernel
#       - gaussian.kernel
#       - laplacian.kernel
#       - cauchy.kernel
#       - inv.multiquadric.kernel
#       - poly.kernel
# p : number of steps of the rw kernel. It is meaningful and is used only with rw.kernel. If 0 (def), it is ignored.
# probs: probabilities relative to the quantiles to be computed (def. 0.1 to 1 by 0.1 steps) to optimally cut the edges. Minimum values of W (usually zeros) are not considered.
#        Note: the numeric vector must be monotonically crescent. If intloo=FALSE, probs represents the quantile that will be used to filter the network.
# opt.fun : function. Function implementing the metric to choice the optimal score threshold for classification. The accuracy (compute.acc) is the default.
#           Available functions:
#                 - compute.F : F-score
#                 - compute.acc : accuracy
#           N.B.: any function having two arguments representing the vector of predicted and true labels can be in principle used.
# score.probs: probabilities relative to the quantiles to be tested (def. 0.1 to 1 by 0.1 steps) to find the optimal score threshold for the classification.
# seed : initialization of the random generator to take into account the randomness due to ties. If seed=0, no initialization (def: 1).
# n.feat : number of selected features based on moderated t-statistic method computed on the training set. 
# ... : optional arguments for the function score. 
# Output: 
# a list with the following elements:
#   - avg.res : a vector with the following components: average AUC, accuracy, F-score, precision, recall across repetitions. 
#   - st.dev : a vector with the standard deviation of the above measures. 
#   - res.full : a matrix with a number of rows equal to the number of repetitions rep, and having as columns AUC, accuracy, F-score, precision, recall obtained at each repetition. 
#   - scores : a list where each component is an integer vector with the scores computed by P-Net on the test set at each repetition of the k-fold CV. 
#   - sel.feat : a list where each component contains k integer vectors with the indices of the selected features at each repetition of the k-fold CV. Indices refer to the rows of M. The length of the list is equal to rep. 
#   - sel.quantile: a list where each component is a vector with k quantiles corresponding to the computed network edge threshold. The length of the list is equal to rep. 
#   - sel.score.quantile : a list where each component is the quantile corresponding to the best score threshold at each repetition. 
#   - pred.round: a list where each component is a vector with the labels predicted on the test set. The length of the list is equal to rep.
setGeneric("pnet.modtstat.repcv", 
                 function(M, ind.pos, kk=5, rep=1, sim=cor, cor.method="pearson", score=eav.score,  ker=list(rw.kernel=rw.kernel), p=0, probs=seq(0.1, 1, 0.1), opt.fun=compute.acc, score.probs=seq(0.01, 1, 0.01), seed=1, n.feat=1000, intloo=TRUE, ...) standardGeneric("pnet.modtstat.repcv"));
				 
setMethod("pnet.modtstat.repcv", signature(M="matrix"), 
  function(M, ind.pos, kk=5, rep=1, sim=cor, cor.method="pearson", score=eav.score,  ker=list(rw.kernel=rw.kernel), p=0, probs=seq(0.1, 1, 0.1), opt.fun=compute.acc, score.probs=seq(0.01, 1, 0.01), seed=1, n.feat=1000, intloo=TRUE, ...) {

    if (seed!=0)
       set.seed(seed);
    n <- nrow(M);    # number of features
    n.ex <- ncol(M); #number of examples 
    x <- 1:n.ex;     

    labels <- integer(n.ex);
    labels[ind.pos] <- 1;

    kmatrices <- vector("list", kk);
    sel.quantile <- vector("list", rep);
    scores <- vector("list", rep);
    thresh <- numeric(kk);
    best.quantile.score <- integer(kk);
    sel.score.quantile <- vector("list", rep);
    sel.feat <- vector("list", rep); 
    pred.round <- vector("list", rep);

    avg.res <- st.dev <- numeric(5);
    names(avg.res) <- names(st.dev) <- c("AUC", "Acc", "F", "Prec", "Rec");
    res.full <- matrix(numeric(length(avg.res)*rep), nrow=rep);
    colnames(res.full) <- names(avg.res);
    score.name <- return.name(score);

  #repeated 5-fold cv
  for (z in 1:rep) {
	f <- do.stratified.cv.data(x, ind.pos, k=kk, seed=seed+z);  

        s  <- numeric(n.ex); 
	names(s) <- colnames(M); 

        # feature selection with Moderated T-test on the training data:
        sorted.index <- vector("list", kk);  
        for (i in 1:kk)  {
	  ind.train <- setdiff(x, c(f$fold.non.positives[[i]],f$fold.positives[[i]]));
          labels.train <- labels[ind.train];
          
          design <- as.factor(labels.train);
          design <- model.matrix(~design)
          fit <- lmFit (M[,ind.train], design);    
          fit <- eBayes (fit);
          p.value <- topTable(fit, adjust = "fdr", number= n.feat);  
          
          sorted.index[[i]] <- match(rownames(p.value), rownames(M));     
        }
        	
        sel.feat[[z]] <- sorted.index;

        # computing the kernel on the selected features (over all the vectors of sequence - 1:kk):
        ker.name <- names(ker);
        for (i in 1:kk)  {
	    if (ker.name == "rw.kernel"  || ker.name == "identity.kernel") {    
               C <- sim(M[sel.feat[[z]][[i]],], method=cor.method);  					    
	       C[C<0 | is.na(C)] <- 0;								    
	       K <- ker[[1]](C); 											    
            if (p>1)  													    
    	       K <- p.step.rw.kernel(K,p); 								    
            }  else 														    
               K <- ker[[1]](t(M[sel.feat[[z]][[i]],]));

        kmatrices[[i]] <- K;
        }

	# computing the optimal threshold for each fold of M (by internal loo)
        q <- numeric(kk);
        if (intloo) {  
	   for (i in 1:kk)  {
	       ind.train <- setdiff(x, c(f$fold.non.positives[[i]],f$fold.positives[[i]]));
	       ind.pos.train<-setdiff(ind.pos,f$fold.positives[[i]]);

               if (score.name == "KNN.score"  || score.name == "NN.score" || score.name == "eav.score" || score.name == "WSLD.score") {
	           q[i] <- optimize.thresh.by.loo(kmatrices[[i]], ind.pos.train, ind.train, score=score, probs=probs, norm=FALSE, ...)$quantile;
               } else
                   q[i] <- optimize.thresh.by.loo(kmatrices[[i]], ind.pos.train, ind.train, score=score, probs=probs, ...)$quantile; 

               thresh[i] <- quantile(kmatrices[[i]], probs=q[i]);
	   }
	} else {
          for (i in 1:kk) {
              thresh[i] <- quantile(kmatrices[[i]],probs=probs);
              q[i]=probs;       
          }	
	}
        
        sel.quantile[[z]] <- q;

        # net filtering for test set
        for (i in 1:kk)  {
            kmatrices[[i]][kmatrices[[i]]< thresh[i]] <- 0;
	}

	# external cross validation (computed on the test set):
	for (i in 1:kk)  {
	  ind.test <- c(f$fold.non.positives[[i]],f$fold.positives[[i]]);
	  ind.pos.train<-setdiff(ind.pos,f$fold.positives[[i]]);

          if (score.name == "KNN.score"  || score.name == "NN.score" || score.name == "eav.score" || score.name == "WSLD.score") {
	      s[ind.test] <- score(kmatrices[[i]], ind.test, ind.pos.train, norm=FALSE, ...);
          } else
              s[ind.test] <- score(kmatrices[[i]], ind.test, ind.pos.train, ...);

	}
	mm<-max(s);
	if (mm > 0) {
	  scores[[z]] <- s / mm;
        }

   ## Estimation of the scores on the training set and of the optimal score theshold by internal loo on the training set
        for (i in 1:kk)  {
	    ind.train <- setdiff(x, c(f$fold.non.positives[[i]],f$fold.positives[[i]]));
	    ind.pos.train<-setdiff(ind.pos,f$fold.positives[[i]]);

            if (score.name == "KNN.score"  || score.name == "NN.score" || score.name == "eav.score" || score.name == "WSLD.score") {
	        score.train<- score(kmatrices[[i]], ind.train, ind.pos.train, norm=FALSE, ...);
            } else
                score.train<- score(kmatrices[[i]], ind.train, ind.pos.train, ...);

            labels.train <- labels[ind.train];
            best.quantile.score[i] <- find.best.quantile.score(score.train, labels.train, opt.fun=opt.fun, score.probs=score.probs);
	}

   # Computing performances on the test set
      AUC <- numeric(kk);
      ACC <- numeric (kk);
      res.partial <- matrix (, nrow=kk, ncol=3); 

      pred.folds <- numeric(kk); 
      for (i in 1:kk) {
          ind.test <- c(f$fold.non.positives[[i]],f$fold.positives[[i]]);
          labels.test <- labels[ind.test];
          AUC[i] <- AUC.single (scores[[z]][ind.test], labels.test);

          sel.score.quantile[[z]] <- mean(best.quantile.score[i]);
          thresh <- quantile(scores[[z]][ind.test], probs=sel.score.quantile[[z]]);
          pred.labels <- ifelse(scores[[z]][ind.test]>=thresh, 1, 0);
          ACC[i] <- compute.acc(pred.labels, labels.test);
          res.partial[i,]<- compute.F.prec.rec(pred.labels, labels.test);

          pred.folds[ind.test] <- pred.labels; 
      }
   
   pred.round[[z]] <- pred.folds;

   res.full[z, "AUC"] <- mean(AUC);
   res.full[z, "Acc"] <- mean(ACC);
   res.full[z, 3:5] <- apply(res.partial ,2,mean);
   cat("End repetition ... ", z, "\n");
   }#end multiple repetitions

   avg.res <- apply(res.full,2,mean);
   st.dev <- apply(res.full,2,sd);
   return (list (avg.res=avg.res, st.dev=st.dev, res.full=res.full, scores=scores, sel.feat=sel.feat, sel.quantile=sel.quantile, sel.score.quantile=sel.score.quantile, pred.round=pred.round)); 

})


# Function "park_wf": Feature selection, through corrected t-test, and examples ranking through multiple k-fold CV 
# (cross-validation).
#
# This function both selects the features using corrected t-test (default is the Bonferroni correction) and ranks examples 
# using P-Net and the selected features. The selection of the features is performed applying the t-test on each feature and
# using the information from all the labeled patients. The selection of the features 
# and the following construction of the Kernel matrix are computed only once for all the repetitions of the k-fold CV.
# Moreover, we exploit the information coming from the unlabeled patients in the construction of the Kernel matrix. 
# Then, the network is filtered with different thresholds (edge.thresh) and all the following steps are repeated on each of
# these filtered networks.
# During the k-fold CV procedure, the ranking is performed with respect to the test set, while the selection of the optimal 
# score threshold is computed on the training set. The scores computed on the test set and 
# the optimal score threshold are exploited to predict the label of each patient. True and predicted labels are compared to compute 
# AUC, accuracy, F-score, precision, recall and specificity.
# The k-fold CV can be repeated and the aforementioned metrics are stored for each repetition on each filtered network. 
# Input:
# M: matrix of examples data. Rows correspond to features (e.g. expression for a given gene), columns to examples (patients).
# ind.pos: integer vector with the indices of the positive examples. Indices correspond to the columns of M.
# ind.unl: integer vector with the indices of the unlabeled examples. Indices correspond to the columns of M. 
# kk: number of folds for the k-fold CV (def. 5).
# rep: number of repetitions of the k-fold CV procedure (def. 1).
# sim: similarity function to be used between patients (def. cor).
# cor.method: method for the choice of the correlation type (def. "pearson"). Other possibilities are "kendall" or "spearman".
# score: function. It must be a kernel-based score method:
#                 - KNN.score 
#                 - NN.score 
#                 - eav.score  (default)
#                 - WSLD.score
# 		  		  - tot.score
# 		  		  - diff.score
#		  		  - dnorm.score
# k: number of neighbours for the KNN.score (def. 5). Meaningful only if the KNN.score is selected.
# d: parameter for the WSLD.score (def. 2). Meaningful only if the WSLD.score is selected.
# ker: list containing the kernel method. Possible kernels:
#       - rw.kernel (default)
#       - identity.kernel
#       - linear.kernel
#       - gaussian.kernel
#       - laplacian.kernel
#       - cauchy.kernel
#       - inv.multiquadric.kernel
#       - poly.kernel
# p: number of steps of the rw kernel. It is meaningful and it is used only with rw.kernel. If 0 (def) it is ignored.
# edge.thresh: vector of probabilities relative to the quantiles to be computed (def. 0.1 to 0.9 by 0.1 steps) to cut the edges. 
#              Minimum values of W (usually zeros) are not considered.
# opt.fun: function. Function implementing the metric to choice the optimal score threshold for classification.
#           Available functions:
#                 - compute.F  : F-score
#                 - compute.acc  : accuracy (def)
#           N.B.: any function having two arguments representing the vector of predicted and true labels can be in principle used.
# score.probs: probabilities relative to the quantiles to be tested (def. 0.1 to 1 by 0.1 steps) to find the optimal score threshold for the classification.
# seed: initialization of the random generator to take into account the randomness due to ties. If seed=0 no initialization (def. 1).
# p.value_thresh: theshold to select the most meaningful features based on corrected t-test computed using all the labeled examples. 
# adj.method: method to adjust the set of p-values returned by t-test (def. "bonferroni"). Available methods are: "holm", "hochberg", "hommel", 
#            "bonferroni" (default), "BH", "BY", "fdr", "none". 
# ...: optional arguments for the function score.
# Output:
# a list with the following elements:
# 	- avg.res: a list of length(edge.thresh) where each component is a vector of metrics (AUC, Acc, Prec, Rec, Spec, Fmeas) averaged across all the repetition 
#              of the k-fold CV for a specific "edge.thresh" value. 
# 	- st.dev: a list of length(edge.thresh) where each element is a vector of standard devations of the above metrics. 
# 	- res.full: a matrix with a number of rows equal to (rep*edge.thresh) and having as columns "Rep" (the considered repetition
#               for a specific edge.thresh), "edge.thresh" used to cut the edges in the net, "AUC", 
#               "accuracy", "F-score", "precision", "recall", "specificity" obtained at each repetition of the k-fold CV for each net threshold.
# 	- scores: a list where each element stores the computed scores (for the labeled examples) for each specific edge.thresh
#             at each repetition.
# 	- sel.feat: vector containing the indices of the selected features.
# 	- best.quantile.score: a list where each component stores the quantile corresponding to the best score threshold
#                          computed for each fold at each repetition for every specific value of "edge.thresh".
# 	- preds: a list where each element stores the labels predicted on the test set at each repetition 
#            for every specific value of "edge.thresh".

setGeneric("park_wf", 
                 function(M, ind.pos, ind.unl, kk=5, rep=1, sim=cor, cor.method="pearson", score=eav.score, k=5, d=2, ker=list(rw.kernel=rw.kernel), p=0, edge.thresh=seq(0.1, 0.9, 0.1), opt.fun=compute.acc, score.probs=seq(0.01, 1, 0.01), seed=1, p.value_thresh=0.01, adj.method="bonferroni", ...) standardGeneric("park_wf"));
				 
setMethod("park_wf", signature(M="matrix"), 
  function(M, ind.pos, ind.unl, kk=5, rep=1, sim=cor, cor.method="pearson", score=eav.score, k=5, d=2, ker=list(rw.kernel=rw.kernel), p=0, edge.thresh=seq(0.1, 0.9, 0.1), opt.fun=compute.acc, score.probs=seq(0.01, 1, 0.01), seed=1, p.value_thresh=0.01, adj.method="bonferroni", ...) {

 if (seed!=0)
 	set.seed(seed);
 n <- nrow(M);     # number of features
 n.ex <- ncol(M);  # number of observations 
 ind.pos <- sort(ind.pos);
 ind.unl <- sort(ind.unl);
 x <- 1:n.ex;      # indices of all the patients
 y <- x[-ind.unl]; # indices of all labeled patients

 labels <- integer(n.ex);
 labels[ind.pos] <- 1;
 labels[ind.unl] <- -1; #-1 stands for unlabeled patient

 score.name <- return.name(score);
 ker.name <- names(ker);
 
 #definition of useful variables to store results:
 res.full <- matrix(, nrow=(rep*length(edge.thresh)), ncol=8); 
 colnames(res.full) <- c("Rep", "edge.thresh", "AUC", "Acc", "Prec", "Rec", "Spec", "Fmeas");
 rownames(res.full) <- 1:(rep*length(edge.thresh));
 scores <- rep(list(rep(list(rep(NA, length(x))), rep)), length(edge.thresh)); 
 sel.score.quantile <- rep(list(list()), rep); 
 preds <- rep(list(list()), length(edge.thresh)); 
 count <- 1;
 
 #feature selection with corrected t-test on all the labeled patients:
 p.value <- integer(n);

 for (j in 1:n) {
        	p.value[j] <- t.test(M[j,y]~labels[-ind.unl], mu=0, alternative="two.sided", conf.level=0.95, var.eq=F, paired=F)$p.value;
 }

 #Correction of the p-values:
 p.value.adj <- p.adjust(p.value, method=adj.method);  
 sorted <- sort.int(p.value.adj, index.return=TRUE); 
 sf <- sorted$ix[which(sorted$x < p.value_thresh)];

 #Computing the kernel on the selected features and all patients (labeled and unlabeled):
 if (ker.name == "rw.kernel" || ker.name == "identity.kernel") {    
 	C <- sim(M[sf,], method=cor.method);  					    
 	C[C<0 | is.na(C)] <- 0;								    
 	K <- ker[[1]](C, ...); 											    
 	if (p>1)  													    
    	K <- p.step.rw.kernel(K,p); 								    
 } else 														    
        K <- ker[[1]](t(M[sf,]), ...);
 
 #Net filtering based on the arbitrary threshold "edge.thresh":
 best.quantile.score <- vector("list", length(edge.thresh)); 

 for (z in 1:length(edge.thresh)) {
    quantiles <- rep(list(rep(NA, kk)), rep); 

 	K1 <- K; 
 	thresh <- quantile(K1, probs=edge.thresh[z]);
 	K1[K1 < thresh] <- 0; 

 	#repeated k-fold cv
 	for (i in 1:rep) {
 		f <- do.stratified.cv.data(y, positives=ind.pos, k=kk, seed=seed+i);  
		
 		#repetition on folds
 		for (w in 1:kk) {

			ind.test <- sort(c(f$fold.non.positives[[w]], f$fold.positives[[w]]));
			ind.train <- setdiff(y, ind.test);
 			ind.pos.train <- setdiff(ind.pos, f$fold.positives[[w]]);
    		labels.train <- labels[ind.train];
			labels.test <- labels[ind.test];
 			
 			#Estimation of the scores on the training set and of the optimal score threshold by
		  	#internal loo on the training set:
 			diag(K1) <- 0; 
			if (score.name == "KNN.score") { 				  
 	    		score.train <- score(K1, x=ind.train, x.pos=ind.pos.train, k=k, norm=FALSE, ...)   
      		} else if (score.name == "WSLD.score") {		  
 	    		score.train <- score(K1, x=ind.train, x.pos=ind.pos.train, d=d, norm=FALSE, ...) 
      		} else if (score.name == "tot.score" || score.name == "diff.score" || score.name == "dnorm.score") {
          		score.train <- score(K1, x=ind.train, x.pos=ind.pos.train) 
      		} else										  
 	    		score.train <- score(K1, x=ind.train, x.pos=ind.pos.train, norm=FALSE, ...);

			quantiles[[i]][w] <- find.best.quantile.score(score.train, labels.train, opt.fun=opt.fun, score.probs=score.probs); 

 			#External cross validation (computed on the test set):
 			if (score.name == "KNN.score") {  				  
 	  	  		scores[[z]][[i]][ind.test] <- score(K1, ind.test, ind.pos.train, k=k, norm=FALSE, ...)
 			} else if (score.name == "WSLD.score") {		  
 				scores[[z]][[i]][ind.test] <- score(K1, ind.test, ind.pos.train, d=d, norm=FALSE, ...) 
 			} else if (score.name == "tot.score" || score.name == "diff.score" || score.name == "dnorm.score") {
        		scores[[z]][[i]][ind.test] <- score(K1, ind.test, ind.pos.train)
 			} else										  
 	  			scores[[z]][[i]][ind.test] <- score(K1, ind.test, ind.pos.train, norm=FALSE, ...);

  		} #end repetition on folds 

 		#Computing performances on the test set:
		mm<-max(scores[[z]][[i]], na.rm=TRUE);
		if (mm > 0) {
			scores[[z]][[i]] <- scores[[z]][[i]] / mm;
		}

		AUC <- numeric (kk);
		ACC <- numeric (kk);
		res.partial <- matrix (, nrow=kk, ncol=4);
		pred.folds <- rep(NA, length(y));

 		for (w in 1:kk) {
			ind.test <- c(f$fold.non.positives[[w]], f$fold.positives[[w]]);
			labels.test <- labels[ind.test];
			AUC[w] <- AUC.single (scores[[z]][[i]][ind.test], labels.test);

 			thresh <- quantile(scores[[z]][[i]][ind.test], probs=quantiles[[i]][w]);
 			pred.labels <- ifelse(scores[[z]][[i]][ind.test]>=thresh, 1, 0);
 			ACC[w] <- compute.acc(pred.labels, labels.test);
 			res.partial[w,]<- F.measure.single(pred.labels, labels.test)[1:4];
 			
			pred.folds[ind.test] <- pred.labels;
 		}

 		res.full[count,"Rep"] <- i;
		res.full[count,"edge.thresh"] <- edge.thresh[z];
		res.full[count,"AUC"] <- mean(AUC);
		res.full[count,"Acc"] <- mean(ACC);
 		res.full[count, 5:8] <- apply(res.partial, 2, mean);

		count <- count + 1;

		preds[[z]][[i]] <- pred.folds;
 
 	cat("End repetition ...", i,  " of ", rep, " (Edge to cut the matrix K is ",edge.thresh[z],")", "\n");
 	} #end of all repetitions 
 
 	best.quantile.score[[z]] <- quantiles;
 } #repetition on weights

 avg.res <- vector("list", length(edge.thresh));
 st.dev <- vector("list", length(edge.thresh));

 for (z in 1:length(edge.thresh)) {
 	avg.res[[z]] <- apply(res.full[which(res.full[,2] == edge.thresh[z]),][,2:8], 2, mean);
	st.dev[[z]] <- apply(res.full[which(res.full[,2] == edge.thresh[z]),][,2:8], 2, sd);
	st.dev[[z]][1] <- edge.thresh[z];
 }

 return (list (avg.res=avg.res, st.dev=st.dev, res.full=res.full, scores=scores, sel.feat=sf, best.quantile.score=best.quantile.score, preds=preds)); 

}) 

#################################################################################
### DEFINITION OF THE FUNCTIONS IMPLEMENTING RANKING AND MARKERSET ASSESSMENT ###
#################################################################################

# Function that implements a non parametric test to assess the statistical signficance of the overall ranking obtained with P-Net through cross-validation.
# The test is based on multiple permutations of the labels of the examples, and a p.value is estimated by comparing the AUC 
# achieved with the true labels with those obtained by random permutations. The AUC is computed considering only the labeled examples and the random permutation of 
# the labels is performed among the labeled patients, without considering the unlabeled ones. 
# The threshold to filter the net for the external cross-validation is obtained by averaging across the thresholds estimated at each internal loo.
# Note that if intloo=FALSE no internal loo is executed and the net is filtered according to the quantile of the argument probs.
# Input:
# W : a symmetric matrix representing similarity between examples. It can be a kernel or a correlation matrix, or any other matrix representing meaningful similarities between examples.
# ind.pos : integer vector with the indices of the positive examples. Indices correspond to the row (column) of W.
# ind.unl : integer vector with the indices of the unlabeled examples. Indices correspond to the row (column) of W.  
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
ranking.assessment.cv_unlabeled <-function(W, ind.pos, ind.unl, kk=5, score=eav.score, probs=seq(0.1, 1, 0.1), rep=1000, intloo=TRUE, seed=NULL, ...) {
  
  set.seed(seed);
  m <- nrow(W);
  # pnet computed on the true labels
  res <- pnet.cv(W, kk=kk, ind.pos=ind.pos, score=score, probs=probs, intloo=intloo, seed=seed, ...)
  s <- res$s[-ind.unl];
  names(s) <- rownames(W)[-ind.unl]; 
  q <- mean(res$q[-ind.unl]);
  target <- numeric(m);
  target[ind.pos] <- 1;
  target <- target[-ind.unl];
  auc <- AUC.single(s, target);
  
  # random permutations
  n.pos <- length(ind.pos);
  x <- 1:m;
  x <- x[-ind.unl];  #indices of positives and negatives (unlabeled not present)
  rank.p.value <- 0.0;
  for (i in 1:rep) {
    ind.pos.random <- sample(x, n.pos); #permutation only on labeled patients
    target <- rep(0,m);
    target[ind.pos.random] <- 1;
    target <- target[-ind.unl];
    res <- pnet.cv(W, kk=kk, ind.pos=ind.pos.random, score=score, probs=q,  intloo=FALSE, ...);
    auc.random <- AUC.single(res$s[-ind.unl], target);
    if (auc.random >= auc)
      rank.p.value <- rank.p.value + 1;
  }
  return(list(rank.p.value=rank.p.value/rep, auc=auc, scores=s));
}


# Function that implements a non parametric test to assess the statistical significance of a given biomarker set through cross validation.
# The test is based on multiple random resamplings of the biomarker set, and a p.value is estimated by comparing the AUC 
# achieved with the assessed biomarker set with those obtained by random resampling. The AUC is computed considering only the labeled examples. 
# The threshold to filter the net for the external cross validation is obtained by averaging across the thresholds estimated at each internal loo.
# Note that if intloo=FALSE no internal loo is executed and the net is filtered according to the quantile of the argument probs.
# Input:
# C : a matrix representing examples and the associated features. Rows correspond to features, and columns to examples.
# ind.pos : integer vector with the indices of the positive examples. Indices correspond to the columns of C.
# ind.unl : integer vector with the indices of the unlabeled examples. Indices correspond to the columns of C. 
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
markerset.assessment.cv_unlabeled <-function(C, ind.pos, ind.unl, ind.markers, kk=5, score=eav.score, kernel=rw.kernel, a=2, p=1, probs=seq(0.1, 1, 0.1), rep=1000, intloo=TRUE, seed=NULL, ...) {
  
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
  s <- res$s[-ind.unl];   #remove indices of unlabeled patients
  names(s) <- rownames(K)[-ind.unl]; 
  q <- mean(res$q[-ind.unl]);
  target <- numeric(m);
  target[ind.pos] <- 1;
  target <- target[-ind.unl];
  auc <- AUC.single(s, target);
  
  # random resampling of markers
  marker.p.value <- 0.0;
  n.markers <- length(ind.markers);
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
    auc.random <- AUC.single(res$s[-ind.unl], target);
    if (auc.random >= auc)
      marker.p.value <- marker.p.value + 1;
  }
  return(list(marker.p.value=marker.p.value/rep, auc=auc, scores=s));
}


######################################################################################
# Functions to implement the classification with SVM through multiple hold-out (MCCV)
######################################################################################

# Function "svm.ttest.mccv": classification with SVM through multiple hold-out (Monte Carlo cross-validation) applying a pre-selection filter with t-test.
# A pre-selection filter is applied on the training set to select the most important features.
# Then this function trains a classifier on the training set (considering only the selected features) and the obtained model is used to make predictions on the test set (always considering only the selected features). The hold-out can be repeated "n" times (rep=n) and the resulting AUC, accuracy, F-score, precision and recall on the test set are stored.
# Input:
# M : matrix of examples data. Rows correspond to features, columns to examples.
# ind.pos : integer vector with the indices of the positive examples. Indices correspond to the columns of M.
# rep : number of repetitions of the hold-out procedure (def. 1).
# train.ratio : ratio of the training set data w.r.t. the total number of examples (def. 0.7).
# kernel : kernel method (def. linear) for the SVM. Possible kernels:
#       - linear (default)
#       - polynomial
#       - radial
#       - sigmoid
# seed : initialization of the random generator to take into account the randomness due to ties. If seed=0, no initialization (def: 1).
# degree : parameter needed for the kernel of type polynomial (def: 3).
# gamma : parameter needed for all kernels except linear (def: 1/ncol(M)).
# coef0 : parameter needed for kernels of type polynomial and sigmoid (def: 0).
# cost : cost of constraints violation (def: 1) — it is the ‘C’-constant of the regularization term in the Lagrange formulation.
# n.feat : number of selected features based on T-test computed on the training set.
# ...: parameters for the svm.
# Output: 
# a list with the following elements:
# avg.res : a vector with the following components: average AUC, accuracy, F-score, precision, recall across repetitions.
# st.dev : a vector with the standard deviation of the above measures.
# res.full : a matrix with a number of rows equal to the number of repetitions rep, and having as columns AUC, accuracy, F-score, precision, recall obtained at each repetition.
# pred : a list whose elements are the computed predictions at each repetition.
# sel.feat : a list where each component is an integer vector with the indices of the selected features. Indices refer to the rows of M. The length of the list is equal to rep.
svm.ttest.mccv <- function (M, ind.pos, rep=1, train.ratio=0.7, kernel="linear", degree=3, gamma=1/ncol(M), coef0=0, cost=1, seed=1, n.feat=1000, ...) {

 library ("e1071");
 if (seed!=0)
   set.seed(seed);
 n <- nrow(M); # number of features
 n.ex <- ncol(M); # number of patients

 x <- 1:n.ex;
 ind.pos <- sort(ind.pos);
 labels <- integer(n.ex);
 labels[ind.pos] <- 1; #Bulding of the labels (0/1)

 avg.res <- st.dev <- numeric(5);
 names(avg.res) <- names(st.dev) <- c("AUC", "Acc", "F", "Prec", "Rec");
 res.full <- matrix(numeric(length(avg.res)*rep), nrow=rep);
 colnames(res.full) <- names(avg.res);
 sel.feat <- vector("list", rep);
 pred <- vector("list", rep);

 # repeated multiple hold-out
 for (i in 1:rep) { 

      f <- do.stratified.hold.out.data(x, ind.pos, train.ratio=train.ratio, seed=seed+i);
      
      #feature selection with T-test on the training data:
      ind.train <- sort(x[c(f$train.pos, f$train.neg)]);
      ind.pos.train<-f$train.pos;
      labels.train <- labels[ind.train];
      p.value <- integer(n);

      for (k in 1:n) #t-test applied only to the training set
             p.value[k] <- t.test( M[k,ind.train]~labels.train, mu=0, alternative="two.sided", conf.level=0.95, var.eq=F, paired=F)$p.value;

      sorted <- sort.int(p.value, index.return=TRUE); 
      sorted.index <- integer(n);
      sf <- sel.feat[[i]] <- sorted.index <- sorted[[2]][1:n.feat];

      #Training of the model and prediction with SVM:
      model <- svm(t(M[sf,ind.train]),labels.train, kernel=kernel, type="C-classification", degree=degree, gamma=gamma, coef0=coef0, ...);

      pred[[i]] <- predict (model, t(M[sf,-ind.train]));

      #Computing performances on the test set
      labels.test <- labels[-ind.train];
      pred1 <- as.numeric(pred[[i]]);
      pred1 <- pred1-1;
      res.full[i, "Acc"] <- compute.acc(pred1, labels.test);
      res.full[i, "AUC"] <- AUC.single(pred1, labels.test); 
      res.full[i, 3:5] <- compute.F.prec.rec(pred1, labels.test);

      cat("End repetition ... ", i, "\n");
 }#End repetitions

 avg.res <- apply(res.full,2,mean);
 st.dev <- apply(res.full,2,sd);

 return(list(avg.res=avg.res, st.dev=st.dev, res.full=res.full, pred=pred, sel.feat=sel.feat));

}


