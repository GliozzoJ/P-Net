# P-Net library

**P-Net** is an R software library implementing the novel algorithm Patient-Net, i.e. Network-based ranking of patients with respect to a given phenotype/outcome. It is a semi-supervised algorithm able to assign to each patient a score related to its odds to show a specific phenotype (e.g. clinical outcome, response to a treatment). The scores can be used to rank the patients and, if we choose a proper ''score threshold'', also for the patients' classification. 

## Brief description of the files contained in the P-Net library

The **P-Net** library contains two files:

1. **_pnet.R_** is the main part of the library that contains all the functions implementing the P-Net algorithm: 

	* *compute.thresholds*: Function to select a set of thresholds for a networks.

	* *optimize.thresh.by.loo*: Method to find the "optimal" filtering threshold by loo.

	* *pnet.loo*: Method that implements the P-Net algorithm using a double loo.

	* *pnet.loo.per.patient.filtering*: Method that implements the P-Net algorithm using a double loo (selection of an optimal threshold for each patient).

	* *pnet.cv*: Method that implements the P-Net algorithm using cross-validation.

	* *pnet.heldout*: Method that implements the P-Net algorithm using held-out.

	* *ranking.assessment.loo*: Function that implements a non parametric test to assess the statistical significance of the overall ranking obtained with P-Net through loo.

	* *ranking.assessment.cv*: Function that implements a non parametric test to assess the statistical significance of the overall ranking obtained with P-Net through cross-validation.

	* *ranking.assessment.heldout*: Function that implements a non parametric test to assess the statistical significance of the overall ranking obtained with P-Net through held out.

	* *markerset.assessment.loo*: Function that implements a non parametric test to assess the statistical significance of a given biomarker set through loo.

	* *markerset.assessment.cv*: Function that implements a non parametric test to assess the statistical significance of a given biomarker set through cross validation.

	* *markerset.assessment.heldout*: Function that implements a non parametric test to assess the statistical significance of a given biomarker set through held out.

	* *pnet.class.heldout*: Method that implements the classifier based on the P-Net algorithm using held-out.

	* *pnet.class.cv*: Method that implements the classifier based on the P-Net algorithm using cross-validation.

	* *pnet.class.loo*: Method that implements the classifier based on the P-Net algorithm using loo.

	* *plot.pnet.graph*: Function to plot results obtained through the P-Net algorithm.

	* *return.name*: Function that returns the string name of the function included in a variable.

	* *find.best.quantile.score*: Function to compute the best score threshold.

	* *do.stratified.hold.out.data*: Function to generate data for the stratified hold.out. 

	* *compute.F.prec.rec*: Function to compute the F-measure, precision and recall for a single class.

	* *tot.score*: Method to compute the Total score for a set of vertices.

	* *diff.score*: Method to compute the Differential score for a set of vertices.

	* *dnorm.score*: Method to compute the Differential Normalized score for a set of vertices.


1. **_experiments.R_** contains a series of functions and methods that implement the workflows used in the paper "Network modeling of patients’ biomolecular profiles for clinical phenotype/outcome prediction" submitted to Scientific Reports - Nature. This is the complete list of the considered workflows:

	* *pnet.ttest.mccv*: computes the feature selection by t-test and the examples ranking using P-Net. The generalization performances of the algorithm are assessed through multiple hold-out (Monte Carlo cross-validation). 

	* *pnet.modtstat.repcv*: computes the feature selection by Moderated t-statistic method and the examples ranking using P-Net. The generalization performances of the algorithm are assessed through multiple k-fold cross validation technique. 

	* *park_wf*: computes the feature selection through corrected t-test and the examples ranking using P-Net. The generalization performances of the algorithm are assessed through multiple k-fold cross validation technique.

	* *ranking.assessment.cv_unlabeled*: Function that implements a non parametric test to assess the statistical signficance of the overall ranking obtained with P-Net through cross-validation. This function takes into account the presence of "unlabeled" patients in the dataset.

	* *markerset.assessment.cv_unlabeled*: Function that implements a non parametric test to assess the statistical significance of a given biomarker set through cross validation. This function takes into account the presence of "unlabeled" patients in the dataset.

	* *svm.ttest.mccv*: computes the classification with SVM (Support Vector Machine) through multiple hold-out (Monte Carlo cross-validation) applying a pre-selection filter with t-test.


**N.B.** A detailed explanation of the algorithms is provided inside the files **_pnet.R_** and **_experiments.R_** as a comment placed above each function/method. 

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE.md](LICENSE.md) file for details.

## Loading

The P-Net library can be easily loaded in the R environment. Just open the R environment in the same folder where you download the **P-Net** library and type:

```
source("pnet.R");
source("experiments.R");

```

The following R packages are required to correctly load the **P-Net** library:

```
library("RANKS");
library("PerfMeas");
library("graph");
library("Rgraphviz");
library("RColorBrewer");
library("e1071");
```

## Example of usage

In this section we provide an example of usage of P-Net on publicly available data.
The first step is to download some pancreatic cancer microarray data from ArrayExpress ([link to data](https://www.ebi.ac.uk/arrayexpress/experiments/E-MEXP-2780/)) and load it in R:

```
library("affy");
load("E-MEXP-2780.eSet.r");
```

Phenotypic and expression data (after normalization by RMA function) are stored in two distinct objects:

```
eset <- rma(study);
M <- exprs(eset);
pdata <- pData(eset);
```

M contains RMA normalized signal intensity data from 30 patients and pData contains experimental and clinical features for each patient.
At this point, we can create two vectors of labels based on the time of survival of the patients. Patients are labelled as _poor prognosis_ (label = 1) if their time of survival is below the median value, otherwise the label is _good prognosis_ (label = 0). 


```
median <- median(pdata$Characteristics.TimeOfSurvival.);
labels <- as.numeric(pdata$Characteristics.TimeOfSurvival. < median);
names(labels) <- colnames(M);
ind.good <- which(labels==0);             
ind.poor <- which(labels==1); 
```

We filter out the genes with low expression and low variance, retaining only the probe sets with the highest mean expression for each gene. This step is not mandatory but we perform it to follow the experimental set-up from the paper "Winter C, Kristiansen G, Kersting S, Roy J, Aust D, Knösel T, et al. (2012) Google Goes Cancer: Improving Outcome Prediction for Cancer Patients by Network-Based Ranking of Marker Genes. PLoS Comput Biol 8(5): e1002511". 

```
mean <- apply(M,1,mean);
M1 <- M[mean>6,];

sd <- apply (M1,1,sd);
M2 <- M1 [sd>0.5,];

library('genefilter')
idx = findLargest(rownames(M2), apply(M2,1,mean), data = "hgu133plus2"); 
fM = M2[idx,];
```

Now we can apply P-Net on the data to predict the labels of our patients.
The first step is to build the **patient similarity matrix**, which provides the similarity between biomolecular profiles of patients exploiting the normalized Pearson correlation (by setting to zero all negative values):

```
C <- cor(M, method="pearson");  			
C[C<0 | is.na(C)] <- 0;
```

Then we build the **kernel matrix** from the **similarity matrix**, choosing a proper kernel able to capture the topological characteristics of the underlying graph. In this case, we selected the random walk kernel with 3-steps:

```
K <- rw.kernel(C);
K <- p.step.rw.kernel(K, p=3); 								    
```

Finally, we can use the function **pnet.loo** to perform the double leave-one-out (loo) procedure. The _internal loo_ selects the threshold tau while the _external loo_ evaluates the performance of the method:

```
pnet_res_loo <- pnet.loo(K, ind.pos=ind.poor, score=eav.score, probs=seq(0.01, 1, 0.01), intloo=TRUE); 
```

The output is a list with the following elements:
- s : a vector with scores computed by loo.
- q : the quantile vector corresponding to the thresholds selected for each patient.

Using the functions **pnet.cv** and **pnet.heldout** we can perform the same task but using an external cross-validation or an held-out to evaluate the performances of the method.

Now we can simply plot the graph coming from the similarity matrix to gain some insight about out patients' network:

```
plot.pnet.graph(C, ind.poor, pnet_res_loo$s, colour="green", file.name="graph.eps", mode="dot", fontsize=12, magn=25)
```
![alt text](https://github.com/GliozzoJ/P-Net/blob/master/graph.eps)

In this plot patients are represented as nodes and similarities between patients through edges. Positive patients, i.e. patients having the phenotype of interest are represented as squares, while the other patients are represented with circles. The scores obtained by the algorithm are represented through shaded colours. Patients with the highest scores are represented through more intense colours, while patients having the lowest scores with light shades of color.

_____________________________________________________________________________________________

P-Net library includes also some functions that implement the experimental workflows employed in the paper associated to this library. For instance, the function **pnet.ttest.mccv** performs a Monte Carlo cross-validation (MCCV) using P-Net as learning method and the t-test as feature selection method. The patients' ranking is performed with respect to the test set, while feature selection and net-threshold selection are performed through efficient internal loo on the training set. In the following line of code we run **pnet.ttest.mccv** to perform 5 repetitions of the MCCV using 70% of patients as training set and selecting the 1000 features with lower p-value:

```
pnet_mccv_res <- pnet.ttest.mccv(fM, ind.pos=ind.poor, train.ratio=0.7, rep=5, sim=cor, cor.method="pearson", score=eav.score, ker=list(rw.kernel=rw.kernel), p=8, probs=seq(0.1,1,0.1), opt.fun=compute.acc, score.probs=seq(0.01, 1, 0.01), seed=123, n.feat=1000);
```

The vector "pnet_mccv_res$avg.res" contains the average AUC, accuracy, F-score, precision and recall across repetitions. The AUC in this case is 0.74.

Another function is **pnet.modtstat.repcv** which performs multiple k-fold cross-validation using P-Net as learning method and the moderated t-statistic as feature selection method. The ranking is performed with respect to the test set, while feature selection and net-threshold selection are performed through efficient internal loo on the training set. In the following line of code we run **pnet.modtstat.repcv** to perform 5 repetitions of the multiple 5-fold cross-validation selecting the 1000 features with lower p-value:

```
pnet_repcv_res <- pnet.modtstat.repcv(fM, ind.pos=ind.poor, kk=5, rep=5, sim=cor, cor.method="pearson", score=KNN.score, ker=list(rw.kernel=rw.kernel), p=8, probs=seq(0.1, 1, 0.1), opt.fun=compute.acc, score.probs=seq(0.01, 1, 0.01), seed=123, n.feat=1000);
```

Also in this case the vector "pnet_repcv_res$avg.res" contains the average AUC, accuracy, F-score, precision, recall across repetitions. The AUC in this case is 0.77.
