#' @importFrom S4Vectors DataFrame
#' @importFrom stats pnorm predict median
.improvedCV2 <- function(x, is.spike, sf.cell=NULL, sf.spike=NULL, log.prior=NULL, 
    use.spikes=FALSE, nbins=20, top.prop=0.01, max.iter=50)
# Fits a spline to the log-CV2 values and computes a p-value for its deviation.
#
# written by Aaron Lun
# created 9 February 2017
{
    # Figuring out what rows to fit to.
    all.genes <- seq_len(nrow(x))
    if (any(is.na(is.spike))) { 
        use.spikes <- TRUE
        is.spike <- all.genes
    } else {
        is.spike <- .subset_to_index(is.spike, x, byrow=TRUE)
    }
    if (length(is.spike) < 2L) {
        stop("need at least 2 spike-ins for trend fitting")
    }

    # Extracting statistics.
    if (is.null(log.prior)) {
        is.cell <- seq_len(nrow(x))[-is.spike]

        if (is.null(sf.cell)) {
            sf.cell <- 1
        } 
        sf.cell <- rep(sf.cell, length.out=ncol(x))
        sf.cell <- sf.cell/mean(sf.cell)

        if (is.null(sf.spike)) {
            sf.spike <- 1
        }
        sf.spike <- rep(sf.spike, length.out=ncol(x))
        sf.spike <- sf.spike/mean(sf.spike)

        spike.stats <- .Call(cxx_compute_CV2, x, is.spike-1L, sf.spike, NULL)
        cell.stats <- .Call(cxx_compute_CV2, x, is.cell-1L, sf.cell, NULL)

        means <- vars <- numeric(nrow(x))
        means[is.cell] <- cell.stats[[1]]
        vars[is.cell] <- cell.stats[[2]]
        means[is.spike] <- spike.stats[[1]]
        vars[is.spike] <- spike.stats[[2]]
    } else {
        log.prior <- as.numeric(log.prior)
        all.stats <- .Call(cxx_compute_CV2, x, all.genes-1L, NULL, log.prior)
        means <- all.stats[[1]]
        vars <- all.stats[[2]]
    }
    cv2 <- vars/means^2

    # Pulling out spike-in values.
    ok.means <- is.finite(means) & means > 0
    to.use <- ok.means & is.finite(cv2) & cv2 > 0
    to.use[-is.spike] <- FALSE

    # Ignoring maxed CV2 values due to an outlier (caps at the number of cells).
    ignored <- cv2 >= ncol(x) - 1e-8
    to.use[ignored] <- FALSE
    use.means <- means[to.use]
    use.cv2 <- cv2[to.use]

    tech.FUN <- .fit_trend_improved(use.means, use.cv2, nbins=nbins, top.prop=top.prop, max.iter=max.iter)
    tech.cv2 <- tech.FUN(means)
    log.cv2 <- log(cv2)
    log.tech.cv2 <- log(tech.cv2)
    tech.sd <- median(abs(log.cv2[to.use] - log.tech.cv2[to.use])) * 1.4826

    # Compute p-values.
    p <- rep(1, length(ok.means))
    p[ok.means] <- pnorm(log.cv2[ok.means], mean=log.tech.cv2[ok.means], sd=tech.sd, lower.tail=FALSE)
    if (!use.spikes) {
        p[is.spike] <- NA
    }

    DataFrame(mean=means, cv2=cv2, trend=tech.cv2, ratio=cv2/tech.cv2,
        p.value=p, FDR=p.adjust(p, method="BH"), row.names=rownames(x))
}

#' @importFrom stats nls median quantile
.fit_trend_improved <- function(x, y, nbins=20, max.iter=3, top.prop=0.01) {
    y0 <- log(y)
    x0 <- log(x)

    by.interval <- cut(x0, breaks=nbins)
    freq <- as.numeric(table(by.interval)[by.interval])
    weights <- 1/freq

    # Rough estimation of initial parameters.
    B <- median(y0 + log(x), na.rm=TRUE)
    A <- median(y0[x >= quantile(x, 1-top.prop)])

    # Least squares fitting on the log-transformed 'y', plus robustness iterations. 
    for (i in seq_len(max.iter)) {
        pfit <- nls(y0 ~ log(exp(A) + exp(B)/x), start=list(A=A, B=B), weights=weights)
        paramFUN <- function(x) { exp(coef(pfit)["A"]) + exp(coef(pfit)["B"])/x }

        r <- abs(y0 - log(paramFUN(x)))
        r <- r/(median(r, na.rm=TRUE) * 6)
        r <- pmin(r, 1)
        new.weights <- (1 - r^3)^3/freq
        if (max(abs(new.weights - weights)) < 1e-8) { 
            break
        }
        weights <- new.weights
    }

    paramFUN
}

#' @export
setGeneric("improvedCV2", function(x, ...) standardGeneric("improvedCV2"))

#' @export
setMethod("improvedCV2", "ANY", .improvedCV2)

#' @importFrom SummarizedExperiment assay
#' @export
setMethod("improvedCV2", "SingleCellExperiment", 
          function(x, spike.type=NULL, ..., assay.type="logcounts", logged=NULL, normalized=NULL) {

    log.prior <- NULL
    if (!is.null(logged)) {
        if (logged) {
            log.prior <- .get_log_offset(x)
        }
    } else {
        if (assay.type=="logcounts") {
            log.prior <- .get_log_offset(x)
        } else if (assay.type!="counts" && assay.type!="normcounts") {
            stop("cannot determine if values are logged")
        }
    }

    if (is.null(normalized)) {
        normalized <- FALSE
        if (assay.type=="logcounts" || assay.type=="normcounts") {
            normalized <- TRUE
        } else if (assay.type!="counts") {
            stop("cannot determine if values are normalized")
        }
    }
    
    prep <- .prepare_cv2_data(x, spike.type=spike.type)
    .improvedCV2(assay(x, i=assay.type), is.spike=prep$is.spike, 
                 sf.cell=prep$sf.cell, sf.spike=prep$sf.spike, log.prior=log.prior, ...)          
})

