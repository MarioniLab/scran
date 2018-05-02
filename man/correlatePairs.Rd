\name{correlatePairs}
\alias{correlatePairs}
\alias{correlatePairs,ANY-method}
\alias{correlatePairs,SingleCellExperiment-method}
\alias{correlateNull}

\title{Test for significant correlations}
\description{Identify pairs of genes that are significantly correlated based on a modified Spearman's rho.}

\usage{
correlateNull(ncells, iters=1e6, block=NULL, design=NULL) 

\S4method{correlatePairs}{ANY}(x, null.dist=NULL, tol=1e-8, iters=1e6, 
    block=NULL, design=NULL, lower.bound=NULL, use.names=TRUE, 
    subset.row=NULL, pairings=NULL, per.gene=FALSE,  cache.size=100L,
    BPPARAM=SerialParam())

\S4method{correlatePairs}{SingleCellExperiment}(x, ..., use.names=TRUE, subset.row=NULL, per.gene=FALSE,
    lower.bound=NULL, assay.type="logcounts", get.spikes=FALSE)
}

\arguments{
\item{ncells}{An integer scalar indicating the number of cells in the data set.}
\item{iters}{An integer scalar specifying the number of values in the null distribution.}
\item{block}{A factor specifying the blocking level for each cell.}
\item{design}{A numeric design matrix containing uninteresting factors to be ignored.}
\item{x}{
    A numeric matrix-like object of log-normalized expression values, where rows are genes and columns are cells.
    Alternatively, a SingleCellExperiment object containing such a matrix.
}
\item{null.dist}{A numeric vector of rho values under the null hypothesis.}
\item{BPPARAM}{A BiocParallelParam object to use in \code{bplapply} for parallel processing.}
\item{tol}{A numeric scalar indicating the maximum difference under which two expression values are tied.}
\item{use.names}{
    A logical scalar specifying whether the row names of \code{x} should be used in the output.
    Alternatively, a character vector containing the names to use.
}
\item{subset.row}{A logical, integer or character vector indicating the rows of \code{x} to use to compute correlations.}
\item{pairings}{A \code{NULL} value indicating that all pairwise correlations should be computed;
or a list of 2 vectors of genes between which correlations are to be computed;
or a integer/character matrix with 2 columns of specific gene pairs - see below for details.}
\item{per.gene}{A logical scalar specifying whether statistics should be summarized per gene.}
\item{lower.bound}{A numeric scalar specifying the theoretical lower bound of values in \code{x}, only used when \code{residuals=TRUE}.}
\item{cache.size}{An integer scalar specifying the number of cells for which ranked expression values are stored in memory.
Smaller values can be used in machines with less memory, at the cost of processing speed.}
\item{...}{Additional arguments to pass to \code{correlatePairs,ANY-method}.}
\item{assay.type}{A string specifying which assay values to use.}
\item{get.spikes}{A logical specifying whether spike-in transcripts should be used.}
}

\details{
The aim of the \code{correlatePairs} function is to identify significant correlations between all pairs of genes in \code{x}.
This allows prioritization of genes that are driving systematic substructure in the data set.
By definition, such genes should be correlated as they are behaving in the same manner across cells.
In contrast, genes driven by random noise should not exhibit any correlations with other genes.

An approximation of Spearman's rho is used to quantify correlations robustly based on ranks.
To identify correlated gene pairs, the significance of non-zero correlations is assessed using a permutation test.
The null hypothesis is that the (ranking of) normalized expression across cells should be independent between genes.
This allows us to construct a null distribution by randomizing (ranked) expression within each gene.

The \code{correlateNull} function constructs an empirical null distribution for rho computed with \code{ncells} cells.
When \code{design=NULL}, this is done by shuffling the ranks, calculating the rho and repeating until \code{iters} values are obtained.
The p-value for each gene pair is defined as the tail probability of this distribution at the observed correlation (with some adjustment to avoid zero p-values).
Correction for multiple testing is done using the BH method.

% Yeah, we could use a t-distribution for this, but the empirical distribution is probably more robust if you have few cells (or effects, after batch correction).

For \code{correlatePairs}, a pre-computed empirical distribution can be supplied as \code{null.dist} if available.
Otherwise, it will be automatically constructed via \code{correlateNull} with \code{ncells} set to the number of columns in \code{x}.
For \code{correlatePairs,SingleCellExperiment-method}, correlations should be computed for normalized expression values in the specified \code{assay.type}. 

The lower bound of the p-values is determined by the value of \code{iters}.
If the \code{limited} field is \code{TRUE} in the returned dataframe, it may be possible to obtain lower p-values by increasing \code{iters}.
This should be examined for non-significant pairs, in case some correlations are overlooked due to computational limitations.
The function will automatically raise a warning if any genes are limited in their significance at a FDR of 5\%.

If \code{per.gene=TRUE}, results are summarized on a per-gene basis.
For each gene, all of its pairs are identified, and the corresponding p-values are combined using Simes' method.
This tests whether the gene is involved in significant correlations to \emph{any} other gene.
Setting \code{per.gene=TRUE} is useful for identifying correlated genes without regard to what they are correlated with (e.g., during feature selection).
}

\value{
For \code{correlateNull}, a numeric vector of length \code{iters} is returned containing the sorted correlations under the null hypothesis of no correlations.
Arguments to \code{design} and \code{residuals} are stored in the attributes.

For \code{correlatePairs} with \code{per.gene=FALSE}, a DataFrame is returned with one row per gene pair and the following fields:
\describe{
\item{\code{gene1, gene2}:}{
    Character or integer fields specifying the genes in the pair.
    If \code{use.names=FALSE}, integers are returned representing row indices of \code{x}, otherwise gene names are returned.
}
\item{\code{rho}:}{A numeric field containing the approximate Spearman's rho.}
\item{\code{p.value, FDR}:}{Numeric fields containing the permutation p-value and its BH-corrected equivalent.}
\item{\code{limited}:}{A logical scalar indicating whether the p-value is at its lower bound, defined by \code{iters}.}
} 
Rows are sorted by increasing \code{p.value} and, if tied, decreasing absolute size of \code{rho}.
The exception is if \code{subset.row} is a matrix, in which case each row in the dataframe correspond to a row of \code{subset.row}.

For \code{correlatePairs} with \code{per.gene=TRUE}, a dataframe is returned with one row per gene.
For each row, the \code{rho} field contains the correlation with the largest magnitude across all gene pairs involving the corresponding gene.
The \code{p.value} field contains the Simes p-value, and the \code{FDR} field contains the corresponding adjusted p-value.
}

\section{Accounting for uninteresting variation}{
If the experiment has known (and uninteresting) factors of variation, these can be included in \code{design} or \code{block}.
\code{correlatePairs} will then attempt to ensure that these factors do not drive strong correlations between genes.
Examples might be to block on batch effects or cell cycle phase, which may have substantial but uninteresting effects on expression.

The approach used to remove these factors depends on whether \code{design} or \code{block} is used.
If there is only one factor, e.g., for plate or animal of origin, \code{block} should be used.
Each level of the factor is defined as a separate group of cells.
For each pair of genes, correlations are computed within each group, and a weighted mean based on the group size) of the correlations is taken across all groups.
The same strategy is used to generate the null distribution where ranks are computed and shuffled within each group.

For experiments containing multiple factors or covariates, a design matrix should be passed into \code{design}.
The \code{correlatePairs} function will fit a linear model to the (log-normalized) expression values.
The correlation between a pair of genes is then computed from the residuals of the fitted model.
Similarly, to obtain a null distribution of rho values, normally-distributed random errors are simulated in a fitted model based on \code{design};
    the corresponding residuals are generated from these errors; and the correlation between sets of residuals is computed at each iteration.

We recommend using \code{block} wherever possible (and it will take priority if both \code{block} and \code{design} are specified).
While \code{design} can also be used for one-way layouts, this is not ideal as it involves more work/assumptions:
\itemize{
\item It assumes normality, during both linear modelling and generation of the null distribution.
This assumption is largely unavoidable for complex designs, where some quantitative constraints are required to remove nuisance effects.
\code{x} should generally be log-transformed here, whereas this is not required for (but does not hurt) the first group-based approach.
\item Residuals computed from expression values equal to \code{lower.bound} are set to a constant value below all other residuals.
This preserves ties between zeroes and avoids spurious correlations between genes due to arbitrary tie-breaking.
The value of \code{lower.bound} should be equal to log-prior count used during the log-transformation.
It is automatically taken from \code{metadata(x)$log.exprs.offset} if \code{x} is a SingleCellExperiment object.
}
}

\section{Gene selection}{
The \code{pairings} argument specifies the pairs of genes to compute correlations for:
\itemize{
\item By default, correlations will be computed between all pairs of genes with \code{pairings=NULL}.
Genes that occur earlier in \code{x} are labelled as \code{gene1} in the output DataFrame.
Redundant permutations are not reported.
\item If \code{pairings} is a list of two vectors, correlations will be computed between one gene in the first vector and another gene in the second vector.
This improves efficiency if the only correlations of interest are those between two pre-defined sets of genes.
Genes in the first vector are always reported as \code{gene1}.
\item If \code{pairings} is an integer/character matrix of two columns, each row is assumed to specify a gene pair.
Correlations will then be computed for only those gene pairs, and the returned dataframe will \emph{not} be sorted by p-value.
Genes in the first column of the matrix are always reported as \code{gene1}.
}

If \code{subset.row} is not \code{NULL}, only the genes in the selected subset are used to compute correlations.
This will iteract properly with \code{pairings}, such that genes in \code{pairings} and not in \code{subset.row} will be ignored.
With \code{correlatePairs,SingleCellExperiment-method}, rows corresponding to spike-in transcripts are also removed by default with \code{get.spikes=FALSE}.
This avoids picking up strong technical correlations between pairs of spike-in transcripts.

We recommend setting  \code{subset.row} to contain only the subset of genes of interest.
This reduces computational time by only testing correlations of interest.
For example, we could select only HVGs to focus on genes contributing to cell-to-cell heterogeneity (and thus more likely to be involved in driving substructure).
There is no need to account for HVG pre-selection in multiple testing, because rank correlations are unaffected by the variance.

Lowly-expressed genes can also cause problems when \code{design} is non-\code{NULL} and \code{residuals=TRUE}.
Tied counts, and zeroes in particular, may not result in tied residuals after fitting of the linear model. 
Model fitting may break ties in a consistent manner across genes, yielding large correlations between genes with many zero counts.
Focusing on HVGs should mitigate the detection of these uninteresting correlations, as genes dominated by zeroes will usually have low variance.
}

\section{Approximating Spearman's rho with tied values}{
As previously mentioned, an approximate version of Spearman's rho is used.
Specifically, untied ranks are randomly assigned to any tied values.
This means that a common empirical distribution can be used for all gene pairs, rather than having to do new permutations for every pair to account for the different pattern of ties.
Generally, this modification has little effect on the results for expressed genes (and in any case, differences in library size break ties for normalized expression values).
Some correlations may end up being spuriously large, but this should be handled by the error control machinery after multiplicity correction.
}

\author{
Aaron Lun
}

\seealso{
\code{\link{bpparam}},
\code{\link{cor}}
}

\references{
Phipson B and Smyth GK (2010).
Permutation P-values should never be zero: calculating exact P-values when permutations are randomly drawn.
\emph{Stat. Appl. Genet. Mol. Biol.} 9:Article 39.

Simes RJ (1986). 
An improved Bonferroni procedure for multiple tests of significance. 
\emph{Biometrika} 73:751-754.
}

\examples{
set.seed(0)
ncells <- 100
null.dist <- correlateNull(ncells, iters=100000)
exprs <- matrix(rpois(ncells*100, lambda=10), ncol=ncells)
out <- correlatePairs(exprs, null.dist=null.dist)
hist(out$p.value) 
}

\keyword{
correlation
}
