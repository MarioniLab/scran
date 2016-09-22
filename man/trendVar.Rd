\name{trendVar}
\alias{trendVar}
\alias{trendVar,matrix-method}
\alias{trendVar,SCESet-method}

\title{Fit a variance trend}
\description{Fit a mean-dependent trend to the gene-specific variances in single-cell RNA-seq data.}

\usage{
\S4method{trendVar}{matrix}(x, trend=c("loess", "semiloess"), 
    span=0.3, family="symmetric", degree=1, start=list(), 
    design=NULL, subset.row=NULL)

\S4method{trendVar}{SCESet}(x, subset.row=NULL, ..., 
    assay="exprs", use.spikes=TRUE)
}

\arguments{
\item{x}{
    A numeric matrix of normalized expression values, where each column corresponds to a cell and each row corresponds to a spike-in transcript.
    Alternatively, a SCESet object that contains such values.
}
\item{trend}{A string indicating whether the trend should be polynomial or loess-based.}
\item{span, family, degree}{Arguments to pass to \code{\link{loess}}.}
\item{start}{A list containing starting values for \code{\link{nls}}.}
\item{design}{A numeric matrix describing the systematic factors contributing to expression in each cell.}
\item{subset.row}{A logical, integer or character scalar indicating the rows of \code{x} to use.}
\item{...}{Additional arguments to pass to \code{trendVar,matrix-method}.}
\item{assay}{A string specifying which assay values to use, e.g., \code{counts} or \code{exprs}.}
\item{use.spikes}{A logical scalar specifying whether the trend should be fitted to variances for spike-in transcripts or endogenous genes.}
}

\details{
The strategy is to fit an abundance-dependent trend to the variance of the log-normalized expression for the spike-in transcripts, using \code{trendVar}.
For SCESet objects, these expression values can be computed by \code{\link{normalize}} after setting the size factors, e.g., with \code{\link{computeSpikeFactors}}.
Log-transformed values are used as these are more robust to genes/transcripts with strong expression in only one or two outlier cells.

The mean and variance of the log-CPMs is calculated for each spike-in transcript, and a trend is fitted to the variance against the mean for all transcripts.
The fitted value of this trend represents technical variability due to sequencing, drop-outs during capture, etc.
Variance decomposition to biological and technical components for endogenous genes can then be performed later with \code{\link{decomposeVar}}.

The design matrix can be set if there are factors that should be blocked, e.g., batch effects, known (and uninteresting) clusters.
Otherwise, it will default to an all-ones matrix, effectively treating all cells as part of the same group.

% This function is separated from decomposeVar in order to allow for different design matrices/cells to be used for trend fitting and variance estimation.
% For example, if you have multiple groups, you can fit the trend to the spike-ins for all groups, but estimate the variance individually for each group.
}

\section{Trend fitting options}{
By default, a robust loess curve is used for trend fitting via \code{\link{loess}}.
This protects against genes with very large or very small variances.
Some experimentation with \code{span}, \code{degree} or \code{family} may be required to obtain satisfactory results.
The fit is also dependent on the quality of the spike-ins -- the fit will obviously be poor if the coverage of all spike-ins is low.

Alternatively, when \code{trend="semiloess"}, a non-linear curve of the form
\deqn{y = \frac{ax}{x^n + b}}{y = ax/(x^n + b)}
is fitted to the variances against the means using \code{\link{nls}}, and a loess curve is then fitted to the log-ratios of the observed to fitted values.
The parametric curve reduces the sharpness of the trend for easier loess fitting.
Conversely, the parametric form is not exact, so the loess curve models any remaining trends in the residuals.

In general, the semi-loess setting tends to give smoother curves than loess alone.
It is more robust to the uneven distribution of spike-in transcripts across the covariate range.
However, it tends to be susceptible to convergence issues, and may require some fiddling with the start values to converge properly.
By default, the start values are \code{a=5}, \code{n=5} and \code{b=1}, which can be altered as named arguments in \code{start}.
}

\section{Additional notes on row selection}{
Spike-in transcripts can be selected in \code{trendVar,SCESet-method} using the \code{use.spikes} method. 
By default, \code{use.spikes=TRUE} which means that only rows labelled as spike-ins with \code{isSpike(x)} will be used.

When spike-ins are not available, \code{trendVar} can also be applied directly to the counts for endogenous genes by setting \code{use.spikes=FALSE} 
    (or by manually supplying a matrix of normalized expression for endogenous genes, for \code{trendVar,matrix-method}).
This assumes that most genes exhibit technical variation and little biological variation, e.g., in a homogeneous population.

If \code{use.spikes=NA}, every row will be used for trend fitting, regardless of whether it corresponds to a spike-in transcript or endogenous gene.
Users can also directly specify which rows to use with \code{subset.row}.
This will override any setting of \code{use.spikes}.
}

\value{
A named list is returned, containing:
\describe{
\item{\code{mean}:}{A numeric vector of mean log-CPMs for all spike-in transcripts.}
\item{\code{var}:}{A numeric vector of the variances of log-CPMs for all spike-in transcripts.}
\item{\code{trend}:}{A function that returns the fitted value of the trend at any mean log-CPM.}
\item{\code{design}:}{A numeric matrix, containing the design matrix that was used.}
}
}

\seealso{
\code{\link{nls}},
\code{\link{loess}},
\code{\link{decomposeVar}},
\code{\link{computeSpikeFactors}},
\code{\link{computeSumFactors}},
\code{\link{normalize}}
}

\author{
Aaron Lun
}

\examples{
set.seed(100)

nspikes <- ncells <- 200
spike.means <- 2^runif(nspikes, 3, 8)
spike.disp <- 100/spike.means + 0.5
spike.data <- matrix(rnbinom(nspikes*ncells, mu=spike.means, size=1/spike.disp), ncol=ncells)

ngenes <- 10000
cell.means <- 2^runif(ngenes, 2, 10)
cell.disp <- 100/cell.means + 0.5
cell.data <- matrix(rnbinom(ngenes*ncells, mu=cell.means, size=1/cell.disp), ncol=ncells)

combined <- rbind(cell.data, spike.data)
colnames(combined) <- seq_len(ncells)
rownames(combined) <- seq_len(nrow(combined))
y <- newSCESet(countData=combined)
y <- calculateQCMetrics(y, list(Spike=rep(c(FALSE, TRUE), c(ngenes, nspikes))))
isSpike(y) <- "Spike"

# Normalizing.
y <- computeSumFactors(y) 
y <- computeSpikeFactors(y, general.use=FALSE)
y <- normalize(y)

# Fitting a trend to the spike-ins.
fit <- trendVar(y)
plot(fit$mean, fit$var)
x <- sort(fit$mean)
lines(x, fit$trend(x), col="red", lwd=2)

# Fitting a trend to the endogenous genes. 
fit <- trendVar(y, use.spikes=FALSE)
plot(fit$mean, fit$var)
x <- sort(fit$mean)
lines(x, fit$trend(x), col="red", lwd=2)
}

\keyword{variance}

