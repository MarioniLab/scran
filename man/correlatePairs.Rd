\name{correlatePairs}
\alias{correlatePairs}
\alias{correlatePairs,matrix-method}
\alias{correlatePairs,SCESet-method}
\alias{correlateNull}

\title{Test for significant correlations}
\description{Identify pairs of genes that are significantly correlated based on a modified Spearman's rho.}

\usage{
correlateNull(ncells, iters=1e6, design=NULL, residuals=FALSE) 

\S4method{correlatePairs}{matrix}(x, null.dist=NULL, design=NULL, BPPARAM=SerialParam(), 
    tol=1e-8, iters=1e6, residuals=FALSE, use.names=TRUE, 
    subset.row=NULL, per.gene=FALSE, lower.bound=NULL)

\S4method{correlatePairs}{SCESet}(x, ..., use.names=TRUE, subset.row=NULL, per.gene=FALSE,
    lower.bound=NULL, assay="exprs", get.spikes=FALSE)
}

\arguments{
\item{ncells}{An integer scalar indicating the number of cells in the data set.}
\item{iters}{An integer scalar specifying the number of values in the null distribution.}
\item{design}{A numeric design matrix containing uninteresting factors to be ignored.}
\item{residuals}{A logical scalar indicating whether correlations should be calculated from residuals when \code{design!=NULL}.}
\item{x}{
    A numeric matrix of log-normalized expression values, where rows are genes and columns are cells.
    Alternatively, a SCESet object containing such a matrix.
}
\item{null.dist}{A numeric vector of rho values under the null hypothesis.}
\item{BPPARAM}{A BiocParallelParam object to use in \code{bplapply} for parallel processing.}
\item{tol}{A numeric scalar indicating the maximum difference under which two expression values are tied.}
\item{use.names}{
    A logical scalar specifying whether the row names of \code{exprs} should be used in the output.
    Alternatively, a character vector containing the names to use.
}
\item{subset.row}{A logical, integer or character vector indicating the rows of \code{x} to use.
Alternatively, a list of 2 such vectors, or a integer/character matrix with 2 columns - see below for details.}
\item{per.gene}{A logical scalar specifying whether statistics should be summarized per gene.}
\item{lower.bound}{A numeric scalar specifying the theoretical lower bound of values in \code{x}, only used when \code{residuals=TRUE}.}
\item{...}{Additional arguments to pass to \code{correlatePairs,matrix-method}.}
\item{assay}{A string specifying which assay values to use.}
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
Otherwise, it will be automatically constructed via \code{correlateNull} with \code{ncells} set to the number of columns in \code{exprs}.
If the \code{limited} field is \code{TRUE} in the returned dataframe, it may be possible to obtain lower p-values by increasing \code{iters}.
This should be examined for non-significant pairs, in case some correlations are overlooked due to computational limitations.
The function will automatically raise a warning if any genes are limited in their significance at a FDR of 5\%.

For \code{correlatePairs,SCESet-method}, correlations should be computed for normalized expression values in the specified \code{assay}. 
By default, rows corresponding to spike-in transcripts are removed with \code{get.spikes=FALSE}. 
This avoids picking up strong technical correlations between pairs of spike-in transcripts.
Users can also set \code{subset.row} to specify which genes to test, which will override any setting of \code{get.spikes}.

If \code{per.gene=TRUE}, results are summarized on a per-gene basis.
For each gene, all of its pairs are identified, and the corresponding p-values are combined using Simes' method.
This tests whether the gene is involved in significant correlations to \emph{any} other gene.
Setting \code{per.gene=TRUE} is useful for identifying correlated genes without regard to what they are correlated with (e.g., during feature selection).
}

\value{
For \code{correlateNull}, a numeric vector of length \code{iters} is returned containing the sorted correlations under the null hypothesis of no correlations.
Arguments to \code{design} and \code{residuals} are stored in the attributes.

For \code{correlatePairs} with \code{per.gene=FALSE}, a dataframe is returned with one row per gene pair and the following fields:
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
No sorting is performed, so if \code{subset.row=NULL}, each row of the returned dataframe corresponds to a row of \code{x}.
For the \code{SCESet} method, if \code{get.spikes=FALSE}, rows corresponding to spike-in features will have \code{NA} entries for most fields.
}

\section{Accounting for uninteresting variation}{
If the experiment has known (and uninteresting) factors of variation, these can be included in \code{design}.
These factors will be regressed out to ensure that they do not drive strong correlations between genes.
Examples might be to block on batch effects or cell cycle phase, which may have substantial but uninteresting effects on expression.

The approach used to remove these factors depends on the design matrix.
If there is only one factor in \code{design}, the levels of the factor are defined as separate groups.
For each pair of genes, correlations are computed within each group, and a weighted mean (based on the group size) of the correlations is taken across all groups.
The same strategy is used to generate the null distribution where ranks are computed and shuffled within each group.

For designs containing multiple factors or covariates, a linear model is fitted to the (log-normalized) expression values with \code{design}.
The correlation between a pair of genes is then computed from the residuals of the fitted model.
Similarly, to obtain a null distribution of rho values, normally-distributed random errors are simulated in a fitted model based on \code{design};
    the corresponding residuals are generated from these errors; and the correlation between sets of residuals is computed at each iteration.

The second model-based approach can also be used for one-way layouts by setting \code{residuals=TRUE}.
By default, this is not turned on as the second approach involves more work/assumptions:
\itemize{
\item It assumes normality, during both linear modelling and generation of the null distribution.
This assumption is largely unavoidable for complex designs, where some quantitative constraints are required to remove nuisance effects.
\code{x} should generally be log-transformed here, whereas this is not required for (but does not hurt) the first group-based approach.
\item Residuals computed from expression values equal to \code{lower.bound} are set to a constant value below all other residuals.
This preserves ties between zeroes and avoids spurious correlations between genes due to arbitrary tie-breaking.
The value of \code{lower.bound} should be equal to log-prior count used during the log-transformation.
It is automatically determined from the \code{logExprsOffset} slot if \code{x} is a SCESet object.
}
}

\section{Gene selection}{
We recommend setting  \code{subset.row} to contain only the subset of highly variable genes (HVGs).
This will focus on genes contributing to cell-to-cell heterogeneity (and thus more likely to be involved in driving substructure).
There is no need to account for HVG pre-selection in multiple testing, because rank correlations are unaffected by the variance.
For more genes, set \code{BPPARAM} to use more workers and reduce computational time.

Lowly-expressed genes can also cause problems when \code{design} is non-\code{NULL} and \code{residuals=TRUE}.
Tied counts, and zeroes in particular, may not result in tied residuals after fitting of the linear model. 
Model fitting may break ties in a consistent manner across genes, yielding large correlations between genes with many zero counts.
Focusing on HVGs should mitigate the detection of these uninteresting correlations, as genes dominated by zeroes will usually have low variance.

If \code{subset.row} is a list of two vectors, correlations will be computed between one gene in the first vector and another gene in the second vector.
This improves efficiency if the only correlations of interest are those between two pre-defined sets of genes.
Alternatively, if \code{subset.row} is an integer/character matrix of two columns, each row is assumed to specify a gene pair.
Correlations will then be computed for only those gene pairs, and the returned dataframe will \emph{not} be sorted by p-value.
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
