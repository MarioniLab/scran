#' @importFrom S4Vectors DataFrame
#' @importFrom stats pnorm predict median
.improvedCV2 <- function(x, is.spike, sf.cell=NULL, sf.spike=NULL, log.prior=NULL, 
    use.spikes=FALSE, bw.adjust=1, top.prop=0.01, max.iter=50)
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

    tech.FUN <- .fit_trend_improved(use.means, use.cv2, adjust=bw.adjust, top.prop=top.prop, max.iter=max.iter)
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

#' @importFrom stats nls median quantile coef
.fit_trend_improved <- function(x, y, adjust=1, max.iter=50, top.prop=0.01) {
    y <- log(y)
    w <- .inverse_density_weights(log(x), adjust=adjust)

    # Rough estimation of initial parameters.
    B <- median(y + log(x), na.rm=TRUE)
    A <- median(y[x >= quantile(x, 1-top.prop)])

    fitFUN <- function(X, Y, W) {
        nls(Y ~ log(exp(A) + exp(B)/X), start=list(A=A, B=B), weights=W)
    }

    predFUN <- function(fit) {
        Aest <- exp(coef(fit)["A"])
        Best <- exp(coef(fit)["B"])
        function(x) { Aest  + Best/x }
    }

    logpredFUN <- function(fit) {
        FUN <- predFUN(fit)
        function(x) log(FUN(x))
    }

    fit <- .robustify_fit(x, y, fitFUN, logpredFUN, weights=w, max.iter=max.iter)
    predFUN(fit) 
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

