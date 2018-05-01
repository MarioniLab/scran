\name{scrannews}
\title{scran News}
\encoding{UTF-8}

\section{Version 1.7.28}{\itemize{
\item Modified decomposeVar() to return statistics (but not p-values) for spike-ins when get.spikes=NA.
Added block= argument for mean/variance calculations within each level of a blocking factor, followed by reporting of weighted averages (using Fisher's method for p-values).
Automatically record global statistics in the metadata of the output for use in combineVar(). 
Switched output to a DataFrame object for consistency with other functions.

\item Fixed testVar() to report a p-value of 1 when both the observed and null variances are zero.

\item Allowed passing of arguments to irlba() in denoisePCA() to assist convergence.
Reported low-rank approximations for all genes, regardless of whether they were used in the SVD.
Deprecated design= argument in favour of manual external correction of confounding effects.
Supported use of a vector or DataFrame in technical= instead of a function.

\item Allowed passing of arguments to prcomp_irlba() in buildSNNGraph() to assist convergence.
Allowed passing of arguments to get.knn(), switched default algorithm back to a kd-tree.

\item Added the buildKNNGraph() function to construct a simple k-nearest-neighbours graph.

\item Fixed a number of bugs in mnnCorrect(), migrated code to C++ and parallelized functions.
Added variance shift adjustment, calculation of angles with the biological subspace.

\item Modified trend specification arguments in trendVar() for greater flexibility.
Switched from ns() to robustSmoothSpline() to avoid bugs with unloaded predict.ns().
Added block= argument for mean/variance calculations within each level of a blocking factor.

\item Added option to avoid normalization in the SingleCellExperiment method for improvedCV2().
Switched from ns() to smooth.spline() or robustSmoothSpline() to avoid bugs.

\item Replaced zoo functions with runmed() for calculating the median trend in DM().

\item Added block= argument to correlatePairs() to calculate correlations within each level of a blocking factor.
Deprecated the use of residuals=FALSE for one-way layouts in design=.
Preserve input order of paired genes in the gene1/gene2 output when pairings!=NULL.

\item Added block= argument to overlapExprs() to calculate overlaps within each level of a blocking factor.
Deprecated the use of residuals=FALSE for one-way layouts in design=.
Switched to automatic ranking of genes based on ability to discriminate between groups.
Added rank.type= and direction= arguments to control ranking of genes.

\item Modified combineVar() so that it is aware of the global stats recorded in decomposeVar().
Absence of global statistics in the input DataFrames now results in an error.
Added option to method= to use Stouffer's method with residual d.f.-weighted Z-scores.
Added weighted= argument to allow weighting to be turned off for equal batch representation.

\item Modified the behaviour of min.mean= in computeSumFactors() when clusters!=NULL.
Abundance filtering is now performed within each cluster and for pairs of clusters, rather than globally.

\item Switched to pairwise t-tests in findMarkers(), rather than fitting a global linear model.
Added block= argument for within-block t-tests, the results of which are combined across blocks via Stouffer's method.
Added lfc= argument for testing against a log-fold change threshold.
Added log.p= argument to return log-transformed p-values/FDRs.
Removed empirical Bayes shrinkage as well as the min.mean= argument.

\item Added the makeTechTrend() function for generating a mean-variance trend under Poisson technical noise.

\item Added the multiBlockVar() function for convenient fitting of multiple mean-variance trends per level of a blocking factor. 

\item Added the clusterModularity() function for assessing the cluster-wise modularity after graph-based clustering.

\item Added the parallelPCA() function for performing parallel analysis to choose the number of PCs.

\item Modified convertT() to return raw counts and size factors for CellDataSet output.

\item Deprecated exploreData(), selectorPlot() in favour of iSEE().
}}

\section{Version 1.6.0}{\itemize{
\item Supported parallelization in buildSNNGraph(), overlapExprs() with BPPARAM options.

\item Forced zero-derived residuals to a constant value in correlatePairs(), overlapExprs().

\item Allowed findMarkers() to return IUT p-values, to identify uniquely expressed genes in each cluster.
Added options to specify the direction of the log-fold changes, to focus on upregulated genes in each cluster.

\item Fixed bug in correlatePairs() when per.gene=TRUE and no spike-ins are available. 
Added block.size argument to control caching.

\item Switched all C++ code to use the beachmat API.
Modified several functions to accept ANY matrix-like object, rather than only base matrix objects.

\item quickCluster() with method="igraph" will now merge based on modularity to satisfy min.size requirements.
Added max.size option to restrict the size of the output clusters.

\item Updated the trendVar() interface with parametric, method arguments.
Deprecated the trend="semiloess" option in favour of parametric=TRUE and method="loess".
Modified the NLS equation to guarantee non-negative coefficients of the parametric trend.
Slightly modified the estimation of NLS starting parameters.
Second d.f. of the fitted F-distribution is now reported as df2 in the output.

\item Modified decomposeVar() to automatically use the second d.f. when test="f".

\item Added option in denoisePCA() to return the number of components or the low-rank approximation.
The proportion of variance explained is also stored as an attribute in all return values.

\item Fixed a variety of bugs in mnnCorrect().
}}

\section{Version 1.4.0}{\itemize{
\item 
Switched default BPPARAM to SerialParam() in all functions.

\item
Added run argument to selectorPlot().
Bug fix to avoid adding an empty list.

\item
Added exploreData() function for visualization of scRNA-seq data.

\item
Minor bug fix to DM() when extrapolation is required.

\item
Added check for centred size factors in trendVar(), decomposeVar() methods.
Refactored trendVar() to include automatic start point estimation, location rescaling and df2 estimation.

\item
Moved spike-in specification to the scater package.

\item
Deprecated isSpike<- to avoid confusion over input/output types.

\item 
Generalized sandbag(), cyclone() to work for other classification problems.

\item
Added test="f" option in testVar() to account for additional scatter.

\item
Added per.gene=FALSE option in correlatePairs(), expanded accepted value types for subset.row.
Fixed an integer overflow in correlatePairs().
Also added information on whether the permutation p-value reaches its lower bound.

\item
Added the combineVar() function to combine results from separate decomposeVar() calls.

\item 
Added protection against all-zero rows in technicalCV2().

\item
Added the improvedCV2() function as a more stable alternative to technicalCV2().

\item 
Added the denoisePCA() function to remove technical noise via selection of early principal components.

\item
Removed warning requiring at least twice the max size in computeSumFactors().
Elaborated on the circumstances surrounding negative size factors.
Increased the default number of window sizes to be examined.
Refactored C++ code for increased speed.

\item
Allowed quickCluster() to return a matrix of ranks for use in other clustering methods.
Added method="igraph" option to perform graph-based clustering for large numbers of cells.

\item 
Added the findMarkers() function to automatically identify potential markers for cell clusters.

\item
Added the overlapExprs() function to compute the overlap in expression distributions between groups.

\item
Added the buildSNNGraph() function to build a SNN graph for cells from their expression profiles.

\item
Added the correctMNN() function to perform batch correction based on mutual nearest neighbors.

\item
Streamlined examples when mocking up data sets.
}}

\section{Version 1.2.0}{\itemize{
\item
Transformed correlations to a metric distance in quickCluster().

\item 
Removed normalize() in favour of scater's normalize().

\item 
Switched isSpike()<- to accept a character vector rather than a logical vector, to enforce naming of spike-in sets.
Also added warning code when the specified spike-in sets overlap.

\item
Allowed compute*Factors() functions to directly return the size factors.

\item 
Added selectorPlot() function for interactive plotting.

\item 
Switched to a group-based weighted correlation for one-way layouts in correlatePairs() and correlateNull(),
and to a correlation of residuals for more complex design matrices.

\item
Added phase assignments to the cyclone() output.

\item
Implemented Brennecke et al.'s method in the technicalCV2() function.

\item
Updated convertTo() to store spike-in-specific size factors as offsets.

\item
Moved code and subsetting into C++ to improve memory efficiency.

\item
Switched to loess-based trend fitting as the default in trendVar(), replaced polynomial with semi-loess fitting.

\item
Added significance statistics to output of decomposeVar(), with only the p-values replaced by NAs for spike-ins.

\item
Updated documentation and tests.
}}

\section{Version 1.0.0}{\itemize{
\item
New package scran, for low-level analyses of single-cell RNA sequencing data.
}}

