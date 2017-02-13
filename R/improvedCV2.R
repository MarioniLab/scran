.improvedCV2 <- function(x, is.spike, sf.cell=NULL, sf.spike=NULL, 
                          log.prior=NULL, df=4, use.spikes=FALSE)
# Fits a spline to the log-CV2 values and computes a p-value for its deviation.
#
# written by Aaron Lun
# created 9 February 2017
# last modified 10 February 2017
{
    # Figuring out what rows to fit to.
    all.genes <- seq_len(nrow(x))
    if (any(is.na(is.spike))) { 
        use.spikes <- TRUE
        is.spike <- all.genes
    } else {
        is.spike <- .subset_to_index(is.spike, x, byrow=TRUE)
    }

    # Extracting statistics.
    if (is.null(log.prior)) {
        is.cell <- seq_len(nrow(x))[-is.spike]
        if (is.null(sf.cell)) sf.cell <- 1
        sf.cell <- rep(sf.cell, length.out=ncol(x))
        if (is.null(sf.spike)) sf.spike <- 1
        sf.spike <- rep(sf.spike, length.out=ncol(x))

        spike.stats <- .Call(cxx_compute_CV2, x, is.spike-1L, sf.spike, NULL)
        if (is.character(spike.stats)) stop(spike.stats)
        cell.stats <- .Call(cxx_compute_CV2, x, is.cell-1L, sf.cell, NULL)
        if (is.character(cell.stats)) stop(cell.stats)

        means <- vars <- numeric(nrow(x))
        means[is.cell] <- cell.stats[[1]]
        vars[is.cell] <- cell.stats[[2]]
        means[is.spike] <- spike.stats[[1]]
        vars[is.spike] <- spike.stats[[2]]
    } else {
        log.prior <- as.numeric(log.prior)
        all.stats <- .Call(cxx_compute_CV2, x, all.genes-1L, NULL, log.prior)
        if (is.character(all.stats)) stop(all.stats) 
        means <- all.stats[[1]]
        vars <- all.stats[[2]]
    }

    cv2 <- vars/means^2
    log.means <- log(means)
    log.cv2 <- log(cv2)

    # Pulling out spike-in values.
    ok.means <- is.finite(log.means)
    to.use <- ok.means & is.finite(log.cv2)
    to.use[-is.spike] <- FALSE
    use.log.means <- log.means[to.use]
    use.log.cv2 <- log.cv2[to.use]

    # Fit a spline to the log-variances, compute p-values.
    fit <- lm(use.log.cv2 ~ ns(use.log.means, df=df))
    tech.sd <- sqrt(mean(fit$effects[-seq_len(fit$rank)]^2))
    tech.log.cv2 <- predict(fit, data.frame(use.log.means=log.means[ok.means]))

    p <- rep(1, length(ok.means))
    p[ok.means] <- pnorm(log.cv2[ok.means], mean=tech.log.cv2, sd=tech.sd, lower.tail=FALSE)
    if (!use.spikes) p[is.spike] <- NA

    tech.cv2 <- rep(NA_real_, length(ok.means))    
    tech.cv2[ok.means] <- exp(tech.log.cv2 + tech.sd^2/2) # correcting for variance
    return(data.frame(mean=means, var=vars, cv2=cv2, trend=tech.cv2, 
                      p.value=p, FDR=p.adjust(p, method="BH"), row.names=rownames(x)))
}

setGeneric("improvedCV2", function(x, ...) standardGeneric("improvedCV2"))

setMethod("improvedCV2", "matrix", .improvedCV2)

setMethod("improvedCV2", "SCESet", function(x, spike.type=NULL, ..., assay="exprs", logged=NULL) {
    prep <- .prepare_cv2_data(x, spike.type=spike.type)
    
    log.prior <- NULL
    if (!is.null(logged)) {
        if (logged) log.prior <- x@logExprsOffset
    } else {
        if (assay=="exprs") {
            log.prior <- x@logExprsOffset
        } else if (assay!="counts") {
            stop("cannot determine if values are logged or counts")
        }
    }

    .improvedCV2(assayDataElement(x, assay), is.spike=prep$is.spike, 
                  sf.cell=prep$sf.cell, sf.spike=prep$sf.spike, log.prior=log.prior, ...)          
})

# library(scran); library(MASS); library(splines)
# sce <- readRDS("brain_data.rds")
# x <- exprs(sce)
# is.spike <- isSpike(sce)
# log.prior <- sce@logExprsOffset
