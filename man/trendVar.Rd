\name{trendVar}
\alias{trendVar}
\alias{trendVar,matrix-method}
\alias{trendVar,SCESet-method}

\title{Fit a variance trend}
\description{Fit a mean-dependent trend to the gene-specific variances in single-cell RNA-seq data.}

\usage{
\S4method{trendVar}{matrix}(x, trend=c("loess", "semiloess"), span=0.3, 
    family="symmetric", degree=1, start=NULL, design=NULL, subset.row=NULL)

\S4method{trendVar}{SCESet}(x, subset.row=NULL, ..., assay="exprs", use.spikes=TRUE)
}

\arguments{
\item{x}{
    A numeric matrix of normalized expression values, where each column corresponds to a cell and each row corresponds to a spike-in transcript.
    Alternatively, a SCESet object that contains such values.
}
\item{trend}{A string indicating whether the trend should be fully loess-based, or a mixture of loess and parametric fitting.}
\item{span, family, degree}{Arguments to pass to \code{\link{loess}}.}
\item{start}{A named list of numeric scalars, containing starting values for parametric fitting with \code{\link{nls}}. 
This is automatically generated with the most suitable values if set to \code{NULL}.}
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

The mean and variance of the normalized log-counts is calculated for each spike-in transcript, and a trend is fitted to the variance against the mean for all transcripts.
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
Reasonably good starting values for \code{a}, \code{n} and \code{b} are chosen automatically but these can be set manually with \code{start}.
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

\section{Warning on size factor centring}{ 
If \code{assay="exprs"}, \code{trendVar,SCESet-method} will attempt to determine if the expression values were computed from counts via \code{\link[scater]{normalize}}.
If so, a warning will be issued if the size factors are not centred at unity.
This is because different size factors are typically used for endogenous genes and spike-in transcripts.
If these size factor sets are not centred at the same value, there will be systematic differences in abundance between these features.
This precludes the use of a spike-in fitted trend with abundances for endogenous genes in \code{\link{decomposeVar}}.

For other expression values and in \code{trendVar,matrix-method}, the onus is on the user to ensure that normalization preserves differences in abundance.
In other words, the scaling factors used to normalize each feature should have the same mean.
This ensures that spurious differences in abundance are not introduced by the normalization process.

% Forcibly recomputing the "exprs" values within trendVar() is arguably safer, as normalize() enforces centring.
% However, it's slightly expensive, hidden from the user and restricted to this function.
% Better to request a re-calculation on the user's side, such that the recalculated values are consistently used here and elsewhere.
%
% The attempted check should only fail in pathological cases:
% - false positives when counts are added to "x" after computing "exprs" (fairly unusual to not put them in as raw data first)
% - false positives when normalize() is run with centre_size_factors=FALSE (don't know why you would do that)
% - false negatives when counts are removed from "x" after computing "exprs" (this is pretty stupid)
% - false negatives when size factors are modified by centring but "exprs" is not re-computed (arguably a problem for the entire workflow)
%
% The other options would involve enforcing validity checks on "exprs" within the SCESet object itself.
% However, this is probably too obtrusive for little gain; only this function requires consistency in the absolute values of the size factors.
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
example(computeSpikeFactors) # Using the mocked-up data 'y' from this example.

# Normalizing (gene-based factors for genes, spike-in factors for spike-ins)
y <- computeSumFactors(y) 
y <- computeSpikeFactors(y, general.use=FALSE)
y <- normalize(y)

# Fitting a trend to the spike-ins.
fit <- trendVar(y)
plot(fit$mean, fit$var)
curve(fit$trend(x), col="red", lwd=2, add=TRUE)

# Fitting a trend to the endogenous genes. 
fit.g <- trendVar(y, use.spikes=FALSE)
plot(fit.g$mean, fit.g$var)
curve(fit.g$trend(x), col="red", lwd=2, add=TRUE)
}

\references{
Lun ATL, McCarthy DJ and Marioni JC (2016). 
A step-by-step workflow for low-level analysis of single-cell RNA-seq data with Bioconductor.
\emph{F1000Res.} 5:2122

}

\keyword{variance}

