\name{Denoise with PCA}
\alias{denoisePCA}
\alias{denoisePCA,matrix-method}
\alias{denoisePCA,SCESet-method}

\title{Denoise expression with PCA}
\description{Denoise log-expression data by removing principal components corresponding to technical noise.}

\usage{
\S4method{denoisePCA}{matrix}(x, technical, design=NULL, subset.row=NULL)

\S4method{denoisePCA}{SCESet}(x, ..., subset.row=NULL, assay="exprs", get.spikes=FALSE)
}

\arguments{
\item{x}{A numeric matrix of log-expression values for \code{denoisePCA,matrix-method}, 
or a SCESet object containing such values for \code{denoisePCA,SCESet-method}.}
\item{technical}{A named numeric vector containing the technical variance for each gene in \code{x}.
Alternatively, a function that accepts the average log-expression and returns the technical variance.}
\item{design}{A numeric matrix containing the experimental design.
If \code{NULL}, all cells are assumed to belong to a single group.}
\item{subset.row}{A logical, integer or character vector indicating the rows of \code{x} to use.
All genes are used by default.}
\item{...}{Further arguments to pass to \code{denoisePCA,matrix-method}.}
\item{assay}{A string specifying which assay values to use.}
\item{get.spikes}{A logical specifying whether spike-in transcripts should be used.
This is ignored if \code{subset.row} is specified.}
}

\value{
For \code{denoisePCA,matrix-method}, a numeric matrix is returned containing the selected PCs (columns) for all cells (rows).

For \code{denoisePCA,SCESet-method}, a SCESet object is returned containing the PCs in the \code{reducedDimension} slot.
}

\details{
The aim of this function is to perform a principal components analysis to reduce random technical noise in the data.
Random noise is uncorrelated across genes and should be captured by later PCs, as the variance in the data explained by any single gene is low.
In contrast, biological substructure should be correlated and captured by earlier PCs, as this explains more variance for sets of genes.
The idea is to discard later PCs to remove technical noise and improve the resolution of substructure.

The choice of the number of PCs to discard is based on the estimates of technical variance in \code{technical}.
This either uses the trend function obtained from \code{\link{trendVar}} or the technical components from \code{\link{decomposeVar}}.
(Obviously, the log-expression values used in those functions should be the same as those used here.)
The total technical variance is estimated by summing the values across genes.

PCs are discarded in order of increasing variance, stopping before the discarded variance is greater than the total technical variance. 
This eliminates the last PCs that should contain only uncorrelated noise.
In this manner, it is possible to make a justified choice for the number of PCs to retain.
The retained PCs should now be enriched for biological signal (both variance and substructure).

If \code{design} is specified, the residuals of a linear model fitted to each gene are computed.
Because variances computed from residuals are usually underestimated, the residuals are scaled up so that their variance is equal to the residual variance of the model fit.
This ensures that the sum of variances is not understated, which would lead to more PCs being discarded than appropriate.
}

\author{
Aaron Lun
}

\seealso{
\code{\link{trendVar}},
\code{\link{decomposeVar}}
}

\examples{
# Mocking up some data.
ngenes <- 1000
is.spike <- 1:100
means <- 2^runif(ngenes, 6, 10)
dispersions <- 10/means + 0.2
nsamples <- 50
counts <- matrix(rnbinom(ngenes*nsamples, mu=means, size=1/dispersions), ncol=nsamples)
rownames(counts) <- paste0("Gene", seq_len(ngenes))

# Fitting a trend.
lcounts <- log2(counts + 1)
fit <- trendVar(lcounts, subset.row=is.spike)
dec <- decomposeVar(lcounts, fit)

# Denoising (not including the spike-ins in the PCA;
# spike-ins are automatically removed with the SCESet method). 
pcs <- denoisePCA(lcounts, technical=fit$trend, subset.row=-is.spike)
dim(pcs)

# Gets the same result as above, if you don't have the function.
pcs2 <- denoisePCA(lcounts, technical=setNames(dec$tech, 
    rownames(dec)), subset.row=-is.spike) # same result
dim(pcs2)

# With a design matrix.
design <- model.matrix(~factor(rep(0:1, length.out=nsamples)))
fit3 <- trendVar(lcounts, design=design, subset.row=is.spike)
pcs3 <- denoisePCA(lcounts, technical=fit3$trend, design=design, subset.row=-is.spike)
dim(pcs3)
}

