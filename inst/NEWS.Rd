\name{scranNews}
\title{scran News}
\encoding{UTF-8}

\section{Version 1.28}{\itemize{
\item Added a \code{restricted=} option to \code{quickSubCluster()} to enable subclustering on specific clusters.
}}

\section{Version 1.20.0}{\itemize{
\item All deprecated functions from the previous release are now defunct.

\item Added a \code{simplify=} option to \code{quickSubCluster()} to get the cluster assignments directly.

\item Deprecated \code{combinePValues()} as this is replaced by \code{metapod::combineParallelPValues()}.

\item \code{getClusteredPCs()} now uses \code{bluster::clusterRows()} by default.

\item \code{decideTestsPerLabel()} now automatically detects \code{pval.field=} if not supplied.

\item Added the \code{clusterCells()} wrapper around \pkg{bluster} functionality.

\item Removed the option to pass a matrix in \code{design=} from \code{pseudoBulkDGE()}.

\item Migrated all normalization-related functions (\code{computeSumFactors()}, \code{calculateSumFactors()},
\code{cleanSizeFactors()} and \code{computeSpikeFactors()}) to a better home in \pkg{scuttle}.
Soft-deprecated existing functions.

\item Modified \code{getTopHVGs()} to accept a SingleCellExperiment and compute the DataFrame with \code{modelGeneVar()}.

\item Added \code{fixedPCA()} to compute a PCA with a fixed number of components, a la \code{scater::runPCA()} (but without requiring \pkg{scater}).

\item Modified \code{denoisePCA()} so that it now complains if \code{subset.row=} is not provided.

\item Modified all \code{pairwise*} functions so that the p-value from \code{direction="any"} is derived from the two p-values from the one-sided tests.
This is necessary for correctness with all choices of \code{lfc=} and \code{block=}, at the cost of conservativeness when \code{block=NULL} and \code{lfc} is large.
}}

\section{Version 1.18.0}{\itemize{
\item Deprecated \code{coassignProbs()} as this is replaced by \code{bluster::pairwiseRand()} 

\item Deprecated \code{boostrapCluster()} as this is replaced by \code{bluster::bootstrapStability()}.

\item Deprecated \code{gene.names=} in the various \code{pairwise*} functions as being out of scope.

\item Added the \code{testLinearModel()} function to obtain inferences from a linear model.

\item Modified \code{pseudoBulkDGE()} to use formulas/functions in the \code{design=} argument.
Allow \code{contrast=} to be a character vector to be run through \code{makeContrasts()}.

\item Added the \code{pseudoBulkSpecific()} function to test for semi-label-specific DEGs in pseudo-bulk analyses.

\item Added the \code{summaryMarkerStats()} function to compute some basic summary statistics for marker filtering.

\item Modified \code{row.data=} in \code{findMarkers()} to support list inputs.
Added a \code{add.summary=} option to easily include summary information.

\item Modified \code{combineVar()} and \code{combineCV2()} to support list inputs.

\item Deprecated \code{doubletCells()} as this is replaced by \code{scDblFinder::computeDoubletDensity()}.

\item Deprecated \code{doubletCluster()} as this is replaced by \code{scDblFinder::findDoubletClusters()}.

\item Deprecated \code{doubletRecovery()} as this is replaced by \code{scDblFinder::recoverDoublets()}.

\item Added sparse-optimized variance calculations to \code{modelGeneVar()}, \code{modelGeneCV2()} and related functions, 
which may result in slight changes to the results due to numeric precision.

\item Exported \code{combineBlocks()} to assist combining of block-wise statistics in other packages.

\item Added \code{lowess=} and \code{density.weights=} options to \code{fitTrendVar()} to rescue overfitted curves.

\item Raised an error in \code{denoisePCA()} upon mismatches in the matrix and technical statistics.
}}

\section{Version 1.16.0}{\itemize{
\item Added the quickSubCluster() function for convenient subclustering.

\item Added the bootstrapCluster() function to convenient bootstrapping of cluster stability.

\item Added the coassignProb() function to compute coassignment probabilities of alternative groupings.

\item combineMarkers() and findMarkers() report a summary effect size for each cluster.

\item Added the multiMarkerStats() function to combine statistics from multiple findMarkers() calls.

\item Added the clusterPurity() function to evaluate cluster purity as a quality measure.

\item Added the pseudoBulkDGE() function to easily and safely perform pseudo-bulk DE analyses.
Also added the decideTestsPerLabel() and summarizeTestsPerLabel() utilities.

\item Added the clusterSNNGraph() and clusterKNNGraph() wrapper functions for easier graph-based clustering.
Provided a k-means pre-clustering option to handle large datasets.
}}

\section{Version 1.14.0}{\itemize{
\item Removed deprecated approximate= and pc.approx= arguments.

\item Removed deprecated batch correction functions.

\item Added option to pairwiseTTests() for standardization of log-fold changes.

\item Changed default BSPARAM= to bsparam() in quickCluster(), denoisePCA(), doubletCells() and build*NNGraph().

\item Added the pairwiseBinom() function for pairwise binomial tests of gene expression.

\item Renamed output fields of pairwiseWilcox() to use AUC for less confusion.
Added the lfc= argument to test against a log-fold change.

\item Added the fitTrendVar(), fitTrendCV2(), modelGeneVar(), modelGeneVarWithSpikes(),
modelGeneCV2(), modelCV2WithSpikes(), fitTrendPoisson() and modelGeneVarByPoisson() functions to model variability.

\item Deprecated the trendVar(), technicalCV2(), improvedCV2(), decomposeVar(), trendVar(), testVar(),
makeTechTrend(), multiBlockVar() and multiBlockNorm() functions.

\item Modified combineVar() to not weight by residual d.f. unless specifically instructed.

\item Added the combineCV2() function to combine separate CV2 modelling results.

\item Added the test.type= argument in findMarkers() to switch between pairwise DE tests.
Added the row.data= argument to easily include row metadata in reordered tables.
Deprecated overlapExprs(), which is replaced by type="wilcox" in findMarkers().

\item Added the getTopMarkers() function to easily retrieve marker lists from pairwise DE results.

\item Added the getTopHVGs() function to easily retrieve HVG sets from variance modelling results.

\item In all functions that accept a block= argument,
any level of the blocking factor that cannot yield a result (e.g., due to insufficient degrees of freedom)
will now be completely ignored and not contribute to any statistic.

\item Added the getDenoisedPCs() function for general-purpose PCA-based denoising on non-SingleCellExperiment inputs.
Converted denoisePCA() to a normal function, removed the method for ANY matrix.
Dropped max.rank= default to 50 for greater speed in most cases.

\item Added the calculateSumFactors() function for general-purpose calculation of deconvolution factors on non-SingleCellExperiment inputs.
Converted computeSumFactors() to a normal function, removed the method for ANY input.
Auto-guess min.mean= based on the average library size.

\item Deprecated all special handling of spike-in rows,
which are no longer necessary when spike-ins are stored as alternative experiments.

\item Deprecated general.use= in computeSpikeFactors(), 
which is no longer necessary when spike-ins are stored as alternative experiments.

\item Deprecated parallelPCA(), which has been moved to the PCAtools package.

\item Modified clusterModularity() to return upper-triangular matrices,
fixing a bug where the off-diagonal weights were split into two entries across the diagonal.
Added the as.ratio= argument to return a matrix of log-ratios.
Renamed the get.values= argument to get.weights=.

\item Simplified density calculation in doubletCells() for greater robustness.

\item Added a method="holm-middle" option to combinePValues(), to test if most individual nulls are true.
Added a min.prop= option to control the definition of "most".

\item Added a pval.type="some" option to combineMarkers(), as a compromise between the two other modes.
Added a min.prop= option to tune stringency for pval.type="some" and "any".

\item Added the getClusteredPCs() function to provide a cluster-based heuristic for choosing the number of PCs.

\item Added the neighborsTo*NNGraph() functions to generate (shared) nearest neighbor graphs from pre-computed NN results.

\item Switched to using only the top 10\% of HVGs for the internal PCA in quickCluster().
}}

\section{Version 1.12.0}{\itemize{
\item Added option in quickCluster() to cluster on log-expression profiles.
Modified defaults to use graph-based clustering on log-expression-derived PCs.

\item Modified default choice of ref.clust= in computeSumFactors().
Degrade quietly to library size factors when cluster is too small for all pool sizes.

\item Minor change to cyclone() random number generation for consistency upon parallelization.

\item Added BPPARAM= to correlateNull() for parallelization.
Minor change in random number generation for consistency upon parallelization.

\item Minor change to parallelPCA() random number generation for consistency upon parallelization.

\item Created correlateGenes() function to compute per-gene correlation statistics.

\item Modified correlatePairs() to compute expected rho after all possible tie-breaking permutations.
Deprecated cache.size= as all ranks are now returned as in-memory representations.
Deprecated per.gene= in favour of an external call to correlateGenes().
Deprecated tol= as ties are now directly determined by rowRanks().

\item Switched to BiocSingular for PCA calculations across all functions.
Deprecated approximate= and pc.approx= arguments in favour of BSPARAM=.

\item Deprecated all batch correction functions to prepare for the migration to batchelor.
}}

\section{Version 1.10.0}{\itemize{
\item Removed selectorPlot(), exploreData() functions in favour of iSEE.
    
\item Fixed underflow problem in mnnCorrect() when dealing with the Gaussian kernel.
Dropped the default sigma= in mnnCorrect() for better default performance.

\item Supported parallelized block-wise processing in quickCluster().
Deprecated max.size= in favour of max.cluster.size= in computeSumFactors().
Deprecated get.ranks= in favour of scaledColRanks().

\item Added max.cluster.size= argument to computeSumFactors().
Supported parallelized cluster-wise processing.
Increased all pool sizes to avoid rare failures if number of cells is a multiple of 5.
Minor improvement to how mean filtering is done for rescaling across clusters in computeSumFactors().
Throw errors upon min.mean=NULL, which used to be valid.
Switched positive=TRUE behaviour to use cleanSizeFactors().

\item Added simpleSumFactors() as a simplified alternative to quickCluster() and computeSumFactors().

\item Added the scaledColRanks() function for computing scaled and centred column ranks.

\item Supported parallelized gene-wise processing in trendVar() and decomposeVar().
Support direct use of a factor in design= for efficiency.

\item Added doubletCluster() to detect clusters that consist of doublets of other clusters.

\item Added doubletCells() to detect cells that are doublets of other cells via simulations.

\item Deprecated rand.seed= in buildSNNGraph() in favour of explicit set.seed() call.
Added type= argument for weighting edges based on the number of shared neighbors.

\item Deprecated rand.seed= in buildKNNGraph().

\item Added multiBlockNorm() function for spike-abundance-preserving normalization prior to multi-block variance modelling.

\item Added multiBatchNorm() function for consistent downscaling across batches prior to batch correction.

\item Added cleanSizeFactors() to coerce non-positive size factors to positive values based on number of detected genes.

\item Added the fastMNN() function to provide a faster, more stable alternative for MNN correction.

\item Added BPPARAM= option for parallelized execution in makeTechTrend().
Added approx.npts= option for interpolation-based approximation for many cells.

\item Added pairwiseTTests() for direct calculation of pairwise t-statistics between groups.

\item Added pairwiseWilcox() for direct calculation of pairwise Wilcoxon rank sum tests between groups.

\item Added combineMarkers() to consolidate arbitrary pairwise comparisons into a marker list.

\item Bugfixes to uses of block=, lfc= and design= arguments in findMarkers().
Refactored to use pairwiseTTests() and combineMarkers() internally.
Added BPPARAM= option for parallelized execution.

\item Refactored overlapExprs() to sort by p-value based on pairwiseWilcox() and combineMarkers().
Removed design= argument as it is not compatible with p-value calculations.

\item Bugfixes to the use of Stouffer's Z method in combineVar().

\item Added combinePValues() as a centralized internal function to combine p-values.
}}

\section{Version 1.8.0}{\itemize{
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

