combineVar <- function(..., method=c("z", "fisher", "simes", "berger")) 
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
        } else if (!identical(ref, rownames(x))) {
            stop("gene identities should be the same for all arguments in '...'")
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

    # Combining the one-sided p-values.
    p.combine <- lapply(all.results, FUN="[[", i="p.value")

    method <- match.arg(method)
    if (method=="z") {
        all.z <- lapply(p.combine, FUN=qnorm)
        Z <- .weighted_average_vals(all.z, resid.df)
        p.final <- pnorm(Z)
    } else if (method=="fisher") {
        all.logp <- lapply(p.combine, FUN=log)
        p.final <- pchisq(-2*Reduce("+", all.logp), df=2*length(all.logp), lower.tail=FALSE)
    } else {
        p.final <- .combine_pvalues(do.call(cbind, p.combine), 
                                    pval.type=ifelse(method=="simes", "any", "all"))
    }

    output <- do.call(DataFrame, output)
    output$p.value <- p.final
    output$FDR <- p.adjust(p.final, method="BH")
    rownames(output) <- ref
    return(output)
}

.extract_weightings <- function(tabs, field) {
    output <- vector("list", length(tabs))
    for (x in seq_along(tabs)) { 
        cur.w <- metadata(tabs[[x]])[[field]]
        if (is.null(cur.w)) {
            warning("inputs should come from decomposeVar()")
            cur.w <- 1
        }
        if (!is.null(names(cur.w))) { 
            cur.w <- cur.w[rownames(tabs[[x]])] 
            if (any(is.na(cur.w))) {
                stop("gene names in 'rownames' not in 'metadata(...)$resid.df'")
            }
            cur.w <- unname(cur.w)
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


