\name{improvedCV2}
\alias{improvedCV2}
\alias{improvedCV2,matrix-method}
\alias{improvedCV2,SCESet-method}

\title{Stably model the technical coefficient of variation}
\description{
Model the decompose coefficient of variation as a function of the mean, and determine the significance of highly variable genes.
This is intended to be a more stable version of \code{\link{technicalCV2}}.
}

\usage{
\S4method{improvedCV2}{matrix}(x, is.spike, sf.cell=NULL, sf.spike=NULL, 
    log.prior=NULL, df=4, robust=FALSE, use.spikes=FALSE)

\S4method{improvedCV2}{SCESet}(x, spike.type=NULL, ..., assay="exprs", logged=NULL)
}

\arguments{
\item{x}{
    A numeric matrix of counts or log-expression values, where each column corresponds to a cell and each row corresponds to a spike-in transcript.
    Alternatively, a SCESet object that contains such values.
}
\item{is.spike}{A vector indicating which rows of \code{x} correspond to spike-in transcripts.}
\item{sf.cell}{A numeric vector containing size factors for endogenous genes.}
\item{sf.spike}{A numeric vector containing size factors for spike-in transcripts.}
\item{log.prior}{A numeric scalar specifying the pseudo-count added prior to log-transformation.
If this is set, \code{x} is assumed to contain log-expression values, otherwise it is assumed to contain counts.}
\item{df}{An integer scalar indicating the number of degrees of freedom for the spline fit with \code{\link[splines]{ns}}.}
\item{robust}{A logical scalar indicating whether robust fitting should be performed with \code{\link[MASS]{rlm}}.}
\item{use.spikes}{A logical scalar indicating whether p-values should be returned for spike-in transcripts.}
\item{spike.type}{A character vector containing the names of the spike-in sets to use.}
\item{...}{Additional arguments to pass to \code{improvedCV2,matrix-method}.}
\item{assay}{A string specifying which assay values to use.}
\item{logged}{A logical scalar indicating if \code{assay} contains log-expression values.
This is automatically determined if \code{assay="counts"} or \code{"exprs"}.}
}

\details{
This function will estimate the squared coefficient of variation (CV2) and mean for each spike-in transcript.
Both values are log-transformed and a mean-dependent trend is fitted to the log-CV2 values, using a linear model with a natural spline of degree \code{df}.
The trend is used to obtain the technical contribution to the CV2 for each gene.
The biological contribution is computed by subtracting the technical contribution from the total CV2.

Deviations from the trend are identified by modelling the CV2 estimates for the spike-in transcripts as log-normally distributed around the fitted trend.
This accounts for sampling variance as well as any variability in the true dispersions (e.g., due to transcript-specific amplification biases).
The p-value for each gene is calculated from a one-sided Z-test on the log-CV2, using the fitted value as the mean and the robust scale estimate as the standard deviation.
A Benjamini-Hochberg adjustment is applied to correct for multiple testing.

If \code{log.prior} is specified, \code{x} is assumed to contain log-expression values.
These are converted back to the count scale prior to calculation of the CV2.
Otherwise, \code{x} is assumed to contain raw counts, which need to be normalized with \code{sf.cell} and \code{sf.spike} prior to calculating the CV2.
Note that both sets of size factors are set to 1 by default if their values are not supplied to the function.

For \code{improvedCV2,matrix-method}, the rows corresponding to spike-in transcripts are specified with \code{is.spike}.
These rows will be used for trend fitting, while all other rows are treated as endogenous genes.
By default, p-values are set to \code{NA} for the spike-in transcripts, such that they do not contribute to the multiple testing correction.
This behaviour can be modified with \code{use.spikes=TRUE}, which will return p-values for all features.

For \code{improvedCV2,SCESet-method}, transcripts from spike-in sets named in \code{spike.type} will be used for trend fitting.
If \code{spike.type=NULL}, all spike-in sets listed in \code{x} will be used.
Size factors for endogenous genes are automatically extracted via \code{\link{sizeFactors}}.
Spike-in-specific size factors for \code{spike.type} are extracted from \code{x}, if available; otherwise they are set to the size factors for the endogenous genes.
Note that the spike-in-specific factors must be the same for each set in \code{spike.type}.

Users can also set \code{is.spike} to \code{NA} in \code{improvedCV2,matrix-method}; or \code{spike.type} to \code{NA} in \code{decomposeCV2,SCESet-method}.
In such cases, all rows will be used for trend fitting, and (adjusted) p-values will be reported for all rows.
This should be used in cases where there are no spike-ins.
Here, the assumption is that most endogenous genes do not exhibit high biological variability and thus can be used to model decompose variation. 
}

\value{
A data frame is returned containing one row per row of \code{x} (including both endogenous genes and spike-in transcripts).
Each row contains the following information:
\describe{
\item{\code{mean}:}{A numeric field, containing mean (scaled) counts for all genes and transcripts.}
\item{\code{var}:}{A numeric field, containing the variances for all genes and transcripts.}
\item{\code{cv2}:}{A numeric field, containing CV2 values for all genes and transcripts.}
\item{\code{trend}:}{A numeric field, containing the fitted value of the trend in the CV2 values.
Note that the fitted value is reported for all genes and transcripts, but the trend is only fitted using the transcripts.}
\item{\code{p.value}:}{A numeric field, containing p-values for all endogenous genes (\code{NA} for rows corresponding to spike-in transcripts).}
\item{\code{FDR}:}{A numeric field, containing adjusted p-values for all genes.}
}
}

\seealso{
\code{\link[splines]{ns}},
\code{\link{technicalCV2}}
}

\author{
Aaron Lun
}

\examples{
# Mocking up some data.
ngenes <- 10000
nsamples <- 50
means <- 2^runif(ngenes, 6, 10)
dispersions <- 10/means + 0.2
counts <- matrix(rnbinom(ngenes*nsamples, mu=means, size=1/dispersions), ncol=nsamples)
is.spike <- logical(ngenes)
is.spike[seq_len(500)] <- TRUE

# Running it directly on the counts.
out <- improvedCV2(counts, is.spike)
head(out)
plot(out$mean, out$cv2, log="xy")
points(out$mean, out$trend, col="red", pch=16, cex=0.5)

# Same again with an SCESet.
rownames(counts) <- paste0("X", seq_len(ngenes))
colnames(counts) <- paste0("Y", seq_len(nsamples))
X <- newSCESet(countData=counts)
X <- calculateQCMetrics(X, list(Spikes=is.spike))
setSpike(X) <- "Spikes"

# Dummying up some size factors (for convenience only, use computeSumFactors() instead).
sizeFactors(X) <- 1
X <- computeSpikeFactors(X, general.use=FALSE)

# Running it.
out <- improvedCV2(X, spike.type="Spikes")
head(out)
}

\keyword{variance}


