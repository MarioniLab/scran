#' @importFrom BiocParallel bpstart bpstop
#' @importFrom scuttle .bpNotSharedOrUp .ranksafeQR
#' @importFrom beachmat rowBlockApply
.compute_mean_var <- function(x, block, design, subset.row, block.FUN, residual.FUN, BPPARAM, ...) {
    if (!is.null(subset.row)) {
        x <- x[subset.row,,drop=FALSE]
    }
    if (.bpNotSharedOrUp(BPPARAM)) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM))
    }

    if (is.null(design)) {
        if (!is.null(block)) { 
            if (ncol(x)!=length(block)) {
                stop("length of 'block' should be the same as 'ncol(x)'")
            }
            block <- as.factor(block)
            bnames <- levels(block)
        } else {
            block <- factor(integer(ncol(x)))
            bnames <- NULL
        }

        ncells <- as.integer(table(block))
        resid.df <- ncells - 1L
        if (all(resid.df <= 0L)){ 
            stop("no residual d.f. in any level of 'block' for variance estimation")
        }

        raw.stats <- rowBlockApply(x,
            FUN=block.FUN, 
            block=as.integer(block) - 1L, 
            nblocks=nlevels(block), 
            ..., 
            BPPARAM=BPPARAM)

        means <- do.call(rbind, lapply(raw.stats, "[[", i=1))
        vars <- do.call(rbind, lapply(raw.stats, "[[", i=2))
        colnames(means) <- colnames(vars) <- names(ncells) <- bnames

    } else {
        if (!is.null(block)) {
            stop("cannot specify 'design' with multi-level 'block'")
        }

        ncells <- nrow(design)
        resid.df <- ncells - ncol(design)
        if (resid.df <= 0L) {
            stop("no residual d.f. in 'design' for variance estimation")
        }
        QR <- .ranksafeQR(design)

        # Calculating the residual variance of the fitted linear model.
        raw.stats <- rowBlockApply(x, 
            FUN=residual.FUN, 
            qr=QR$qr, 
            qraux=QR$qraux, 
            ..., 
            BPPARAM=BPPARAM)

        means <- matrix(unlist(lapply(raw.stats, FUN="[[", i=1)))
        vars <- matrix(unlist(lapply(raw.stats, FUN="[[", i=2)))
    }

	rownames(means) <- rownames(vars) <- rownames(x)
    list(means=means, vars=vars, ncells=ncells)
}

dummy.trend.fit <- list(trend=function(x) { rep(NA_real_, length(x)) }, std.dev=NA_real_)

#' @importFrom stats pnorm p.adjust
#' @importFrom S4Vectors DataFrame metadata<-
.decompose_log_exprs <- function(x.means, x.vars, fit.means, fit.vars, ncells, ...) {
    collected <- vector("list", ncol(x.means))
    for (i in seq_along(collected)) {
        fm <- fit.means[,i]
        fv <- fit.vars[,i]
        if (ncells[i] >= 2L) {
            fit <- fitTrendVar(fm, fv, ...)
        } else {
            fit <- dummy.trend.fit
        }

        xm <- unname(x.means[,i])
        xv <- unname(x.vars[,i])
        output <- DataFrame(mean=xm, total=xv, tech=fit$trend(xm))
        output$bio <- output$total - output$tech
        output$p.value <- pnorm(output$bio/output$tech, sd=fit$std.dev, lower.tail=FALSE)
        output$FDR <- p.adjust(output$p.value, method="BH")

        rownames(output) <- rownames(x.means)
        metadata(output) <- c(list(mean=fm, var=fv), fit)
        collected[[i]] <- output
    }
    names(collected) <- colnames(x.means)
    collected
}

#' @importFrom stats pnorm p.adjust
#' @importFrom S4Vectors DataFrame metadata<-
.decompose_cv2 <- function(x.means, x.vars, fit.means, fit.vars, ncells, ...) {
    collected <- vector("list", ncol(x.means))
    for (i in seq_along(collected)) {
        fm <- fit.means[,i]
        fcv2 <- fit.vars[,i]/fm^2
        if (ncells[i] >= 2L) {
            fit <- fitTrendCV2(fm, fcv2, ncells[i], ...)
        } else {
            fit <- dummy.trend.fit
        }

        xm <- unname(x.means[,i])
        xcv2 <- unname(x.vars[,i])/xm^2
        output <- DataFrame(mean=xm, total=xcv2, trend=fit$trend(xm))

        output$ratio <- output$total/output$trend
        output$p.value <- pnorm(output$ratio, mean=1, sd=fit$std.dev, lower.tail=FALSE)
        output$FDR <- p.adjust(output$p.value, method="BH")

        rownames(output) <- rownames(x.means)
        metadata(output) <- c(list(mean=fm, cv2=fcv2), fit)
        collected[[i]] <- output
    }
    names(collected) <- colnames(x.means)
    collected
}

.combine_blocked_statistics <- function(collected, method, equiweight, ncells, geometric=FALSE,
    fields=c("mean", "total", "tech", "bio"), pval="p.value")
{
    combineBlocks(collected, 
        method=method, 
        equiweight=equiweight,
        weights=ncells,
        valid=ncells >= 2L,
        geometric=geometric, 
        ave.fields=fields,
        pval.field=pval)
}

#' @importFrom BiocParallel SerialParam
#' @importFrom scuttle librarySizeFactors
.compute_var_stats_with_spikes <- function(x, spikes, size.factors=NULL, spike.size.factors=NULL, 
    subset.row=NULL, block=NULL, BPPARAM=SerialParam(), ...)
{
    if (is.null(size.factors)) {
        size.factors <- librarySizeFactors(x, subset_row=subset.row)
    } else {
        # Centering just in case; this should almost certainly be true anyway.
        size.factors <- size.factors/mean(size.factors)
    }
    stats.out <- .compute_mean_var(x, block=block, subset.row=subset.row, 
        BPPARAM=BPPARAM, sf=size.factors, ...)

    if (is.null(spike.size.factors)) {
        spike.size.factors <- librarySizeFactors(spikes) # no subset_row here, as that only applies to 'x'.
    }

    # Rescaling so that the mean spike.size.factors is the same as each
    # size.factors in each block.  Note that we do not recenter the
    # size.factors themselves within each block, so as to ensure we are
    # modelling the variance in the same log-expression values that will be
    # used in downstream analyses.
    if (is.null(block)) {
        spike.size.factors <- spike.size.factors / mean(spike.size.factors) # assume mean(size.factors)=1, see above.
    } else {
        by.block <- split(seq_along(block), block)
        for (i in by.block) {
            current <- spike.size.factors[i]
            spike.size.factors[i] <- current / mean(current) * mean(size.factors[i])
        }
    }

    spike.stats <- .compute_mean_var(spikes, block=block, subset.row=subset.row, 
        BPPARAM=BPPARAM, sf=spike.size.factors, ...)

    list(x=stats.out, spikes=spike.stats)
}

#' @importFrom stats density approx
.inverse_density_weights <- function(x, adjust=1) {
    out <- density(x, adjust=adjust, from=min(x), to=max(x))
    w <- 1/approx(out$x, out$y, xout=x)$y 
    w/mean(w)
}

#' @importFrom limma weighted.median
.correct_logged_expectation <- function(x, y, w, FUN) 
# Adjusting for any scale shift due to fitting to the log-values.
# The expectation of the log-values should be the log-expectation
# plus a factor that is dependent on the variance of the raw values
# divided by the squared mean, using a second-order Taylor approximation. 
# If we assume that the standard deviation of the variances is proportional
# to the mean variances with the same constant across all abundances,
# we should be able to correct the discrepancy with a single rescaling factor. 
{
    leftovers <- y/FUN(x)
    med <- weighted.median(leftovers, w, na.rm=TRUE)

    OUT <- function(x) { 
        output <- FUN(x) * med
        names(output) <- names(x)
        output
    }

    # We assume ratios are normally distributed around 1 with some standard deviation.
    std.dev <- unname(weighted.median(abs(leftovers/med - 1), w, na.rm=TRUE)) * 1.4826 
    list(trend=OUT, std.dev=std.dev)
}
