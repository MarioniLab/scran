combineVar <- function(..., method=c("fisher", "simes", "berger")) 
# Combines decomposeVar() results, typically from multiple batches.
#
# written by Aaron Lun
# created 19 January 2017
{
    all.results <- list(...)
    ref <- NULL
    for (x in all.results) {
        if (is.null(ref)) {
            ref <- rownames(x)
        } else {
            stopifnot(identical(ref, rownames(x)))
        }
    }

    # Averaging each value, weighted by number of cells/residual d.f.
    to.average <- c("mean", "total", "bio", "tech")
    output <- vector("list", length(to.average))
    names(output) <- to.average

    all.means <- lapply(all.results, FUN="[[", i="mean")
    n.cells <- .extract_weightings(all.results, "num.cells")
    output$mean <- .weighted_average_vals(all.means, n.cells)

    resid.df <- .extract_weightings(all.results, "resid.df")
    for (val in to.average[-1]) { 
        cur.vals <- lapply(all.results, FUN="[[", i=val)
        output[[val]] <- .weighted_average_vals(cur.vals, resid.df)
    }

    # Combining the p-values.
    p.combine <- lapply(all.results, FUN="[[", i="p.value")
    p.combine <- do.call(cbind, p.combine)

    method <- match.arg(method)
    if (method=="fisher") {
        logp <- -2*rowSums(log(p.combine))
        p.final <- pchisq(logp, df=2*ncol(p.combine), lower.tail=FALSE)
    } else if (method=="simes") {
        p.final <- apply(p.combine, 1, FUN=function(p) { min(p.adjust(p, method="BH")) })
    } else if (method=="berger") {
        p.final <- apply(p.combine, 1, FUN=max)
    }

    output <- do.call(data.frame, output)
    output$p.value <- p.final
    output$FDR <- p.adjust(p.final, method="BH")
    rownames(output) <- ref
    return(output)
}

.extract_weightings <- function(tabs, field) {
    output <- vector("list", length(tabs))
    for (x in seq_along(tabs)) { 
        cur.w <- attr(tabs[[x]], "decomposeVar.stats")[[field]]
        if (is.null(cur.w)) {
            stop("inputs should come from decomposeVar() with store.stats=TRUE")
        }
        output[[x]] <- cur.w 
    } 
    names(output) <- names(tabs)
    return(output)
}

.weighted_average_vals <- function(vals, weights) {
    combined <- total.weight <- 0
    for (x in seq_along(vals)) { 
        cur.weights <- weights[[x]]
        product <- vals[[x]] * cur.weights
        product[cur.weights==0] <- 0 # avoid problems with NA variances.
        combined <- combined + product
        total.weight <- total.weight + cur.weights
    }
    combined/total.weight
}


