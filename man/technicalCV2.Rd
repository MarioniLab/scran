\name{technicalCV2}
\alias{technicalCV2}
\alias{technicalCV2,matrix-method}
\alias{technicalCV2,SCESet-method}

\title{Model the technical coefficient of variation}
\description{Model the technical coefficient of variation as a function of the mean, and determine the significance of highly variable genes.}

\usage{
\S4method{technicalCV2}{matrix}(x, is.spike, sf.cell=NULL, sf.spike=NULL, 
    cv2.limit=0.3, cv2.tol=0.8, min.bio.disp=0.25) 

\S4method{technicalCV2}{SCESet}(x, spike.type=NULL, ..., assay="counts")
}

\arguments{
\item{x}{
    A numeric matrix of counts, where each column corresponds to a cell and each row corresponds to a spike-in transcript.
    Alternatively, a SCESet object that contains such values.
}
\item{is.spike}{A vector indicating which rows of \code{x} correspond to spike-in transcripts.}
\item{sf.cell}{A numeric vector containing size factors for endogenous genes.}
\item{sf.spike}{A numeric vector containing size factors for spike-in transcripts.}
\item{cv2.limit, cv2.tol}{Numeric scalars that determine the minimum mean abundance for the spike-in transcripts to be used for trend fitting.}
\item{min.bio.disp}{A numeric scalar specifying the minimum biological dispersion.}
\item{spike.type}{A character vector containing the names of the spike-in sets to use.}
\item{...}{Additional arguments to pass to \code{technicalCV2,matrix-method}.}
\item{assay}{A string specifying which assay values to use.}
}

\details{
This function will estimate the squared coefficient of variation (CV2) and mean for each spike-in transcript.
A mean-dependent trend is fitted to the CV2 values for the transcripts using a Gamma GLM with \code{\link[statmod]{glmgam.fit}}.
Only high-abundance transcripts are used for stable trend fitting.
(Specifically, a mean threshold is selected by taking all transcripts with CV2 above \code{cv2.limit}, and taking the quantile of this subset at \code{cv2.tol}.
A warning will be thrown and all spike-ins will be used if the subset is empty.)

The trend is used to determine the technical CV2 for each endogenous gene based on its mean.
To identify highly variable genes, the null hypothesis is that the total CV2 for each gene is less than or equal to the technical CV2 plus \code{min.bio.disp}.
Deviations from the null are identified using a chi-squared test.
The additional \code{min.bio.disp} is necessary for a ratio-based test, as otherwise genes with large relative (but small absolute) CV2 would be favoured.

% For NB-distributed counts, log-counts are more normal-looking than the counts themselves.
% (Assessed based on the former yielding variance estimates that are more chi-squared, especially at high dispersions.)
% However, this is probably academic as they're both fairly bad (4-fold for log, 20-fold for row over the type I error rate for a chi-squared test).
% This would be compounded by instability in the curve fitting with limited numbers of spike-ins.

For \code{technicalCV2,matrix-method}, the rows corresponding to spike-in transcripts are specified with \code{is.spike}.
These rows will be used for trend fitting, while all other rows are treated as endogenous genes.
If either \code{sf.cell} or \code{sf.spike} are not specified, the \code{\link[DESeq2]{estimateSizeFactorsForMatrix}} function is applied to compute size factors.

For \code{technicalCV2,SCESet-method}, transcripts from spike-in sets named in \code{spike.type} will be used for trend fitting.
If \code{spike.type=NULL}, all spike-in sets listed in \code{x} will be used.
Size factors for the endogenous genes are automatically extracted via \code{\link{sizeFactors}}.
Spike-in-specific size factors for \code{spike.type} are extracted from \code{x}, if available; otherwise they are set to the size factors for the endogenous genes.
Note that the spike-in-specific factors must be the same for each set in \code{spike.type}.

Users can also set \code{is.spike} to \code{NA} in \code{technicalCV2,matrix-method}; or \code{spike.type} to \code{NA} in \code{technicalCV2,SCESet-method}.
In such cases, all rows will be used for trend fitting, and (adjusted) p-values will be reported for all rows.
This should be used in cases where there are no spike-ins.
Here, the assumption is that most endogenous genes do not exhibit high biological variability and thus can be used to model technical variation. 
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
\code{\link[statmod]{glmgam.fit}},
\code{\link[DESeq2]{estimateSizeFactorsForMatrix}}
}

\author{
Aaron Lun,
based on code from Brennecke et al. (2013)
}

\examples{
# Mocking up some data.
ngenes <- 10000
means <- 2^runif(ngenes, 6, 10)
dispersions <- 10/means + 0.2
nsamples <- 50
counts <- matrix(rnbinom(ngenes*nsamples, mu=means, size=1/dispersions), ncol=nsamples)
is.spike <- logical(ngenes)
is.spike[seq_len(500)] <- TRUE

# Running it directly on the counts.
out <- technicalCV2(counts, is.spike)
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
out <- technicalCV2(X, spike.type="Spikes")
head(out)
}

\keyword{variance}

\references{
Brennecke P, Anders S, Kim JK et al. (2013).
Accounting for technical noise in single-cell RNA-seq experiments.
\emph{Nat. Methods} 10:1093-95
}
