\name{correlatePairs}
\alias{correlatePairs}
\alias{correlatePairs,matrix-method}
\alias{correlatePairs,SCESet-method}
\alias{correlateNull}

\title{Test for significant correlations}
\description{Identify pairs of genes that are significantly correlated based on a modified Spearman's rho.}

\usage{
correlateNull(ncells, iters=1e6, design=NULL) 

\S4method{correlatePairs}{matrix}(x, null.dist=NULL, design=NULL, BPPARAM=bpparam(), 
    use.names=TRUE, tol=1e-8)

\S4method{correlatePairs}{SCESet}(x, ..., assay="exprs", get.spikes=FALSE)
}

\arguments{
\item{ncells}{An integer scalar indicating the number of cells in the data set.}
\item{iters}{An integer scalar specifying the number of values in the null distribution.}
\item{design}{A numeric design matrix describing fixed effects to factorize out.}
\item{x}{
    A numeric matrix of normalized expression values, where rows are genes and columns are cells.
    Alternatively, a SCESet object containing such a matrix.
}
\item{null.dist}{A numeric vector of rho values under the null hypothesis.}
\item{BPPARAM}{A BiocParallelParam object to use in \code{bplapply} for parallel processing.}
\item{use.names}{
    A logical scalar specifying whether the row names of \code{exprs} should be used in the output.
    Alternatively, a character vector containing the names to use.
}
\item{tol}{A numeric scalar indicating the maximum difference under which two expression values are tied.}
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
This is done by shuffling the ranks, calculating the rho and repeating until \code{iters} values are obtained.
The p-value for each gene pair is defined as the tail probability of this distribution at the observed correlation (with some adjustment to avoid zero p-values).
Correction for multiple testing is done using the BH method.

% Yeah, we could use a t-distribution for this, but the empirical distribution is probably more robust if you have few cells (or effects, after batch correction).

For \code{correlatePairs}, a pre-computed empirical distribution can be supplied as \code{null.dist} if available.
Otherwise, it will be automatically constructed via \code{correlateNull} with \code{ncells} set to the number of columns in \code{exprs}.

For \code{correlatePairs,SCESet-method}, correlations should be computed for normalized expression values in the specified \code{assay}. 
By default, rows corresponding to spike-in transcripts are removed with \code{get.spikes=FALSE}. 
This avoids picking up strong technical correlations between pairs of spike-in transcripts.
}

\value{
For \code{correlateNull}, a numeric vector of length \code{iters} is returned containing the sorted correlations under the null hypothesis of no correlations.

For \code{correlatePairs}, a dataframe is returned with one row per gene pair and the following fields:
\describe{
\item{\code{gene1, gene2}:}{
    Character or integer fields specifying the genes in the pair.
    If \code{use.names=FALSE}, integer row indices are returned, otherwise gene names are returned.
}
\item{\code{rho}:}{A numeric field containing the approximate Spearman's rho.}
\item{\code{p.value, FDR}:}{Numeric fields containing the permutation p-value and its BH-corrected equivalent.}
} 
Rows are sorted by increasing \code{p.value} and, if tied, decreasing absolute size of \code{rho}.
}

\section{Gene selection and experimental design}{
Users should select their genes in \code{x} with some care.
Using a top set of 100-200 highly variable genes (HVGs) is recommended.
This will focus on genes contributing to cell-to-cell heterogeneity (and thus more likely to be involved in driving substructure).
There is no need to account for HVG pre-selection in multiple testing, because rank correlations are unaffected by the variance.
For more genes, set \code{BPPARAM} to use more workers and reduce computational time.

If the experiment has known (and uninteresting) factors of variation, e.g., batch effects, cell cycle phase, these can be included in the \code{design}.
Correlations between genes will then be computed using the residual effects of a linear model fitted to the normalized expression values with \code{design}.
Similarly, the null distribution of rho values will be constructed with \code{ncells} set as the residual degrees of freedom in \code{design}.
This procedure ensures that the uninteresting factors do not drive strong correlations between genes.

% Probably okay - main effects in 'x' aside, any scalar/multiplier to 'y' will be preserved in qty, so direction and quality of correlation will be preserved.
% Not exactly equal, though, as we lose a couple of observations when converting to effects.

Note that tied counts may not lead to tied effects, especially for design matrices with real-valued covariates.
As a result, large correlations may be obtained for genes with few non-zero counts.
Some protection is provided by focusing on HVGs, because genes dominated by zeroes will usually have low variance.

% Consider a real-valued covariate of 1:100 for 100 cells, where the count for the last cell is unity in each of two genes, and zeroes everywhere else.
% You'll get correlations of 1 because the residual effects will be negative for the zeros.
% This seems unavoidable; after all, the likelihood of obtaining a zero does change across cells, given the changes in the fitted value you want to correct for.
%
% Don't bother trying to fit a linear model to the ranks, either.
% I thought it would be a generalization of the definition of Spearman's (Pearson's on ranks).
% However, there's no guarantee that unevenly-spaced covariates (or factors, for that matter) will make sense when fitted to ranks.
}

\section{Approximating Spearman's rho with tied values}{
As previously mentioned, an approximate version of Spearman's rho is used.
Specifically, untied ranks are randomly assigned to any tied values.
This means that a common empirical distribution can be used for all gene pairs, rather than having to do new permutations for every pair to account for the different pattern of ties.
Generally, this modification has little effect on the results for expressed genes (and in any case, differences in library size break ties for normalized expression values).
Some correlations may end up being spuriously large, but this should be handled by the error control machinery after multiplicity correction.

For example, counts of zero will have the same normalized log-expression, even if the library sizes are different.
This is because the added prior count is scaled by the library size in \code{\link{cpm.default}}, such that the effect of library size cancels out.
Thus, all zeroes will have tied ranks (with numerical imprecision handled by \code{tol}) and will not inflate the correlations.
For non-zero counts, correlations may be driven by library size differences between cells.
This is, perhaps, less problematic, as two genes with the same counts in both small and large libraries provide evidence for association.
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
