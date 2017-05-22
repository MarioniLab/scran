\name{scrannews}
\title{scran News}
\encoding{UTF-8}

\section{Version 1.5.3}{\itemize{
\item Parallelized buildSNNGraph(), overlapExprs() with BPPARAM options.

\item Forced zero-derived residuals to a constant value in correlatePairs(), overlapExprs().

\item Allowed findMarkers() to return IUT p-values, to identify uniquely expressed genes in each cluster.

\item Fixed bug in correlatePairs() when per.gene=TRUE and no spike-ins are available. 
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

