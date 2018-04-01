#' @importFrom stats pnorm pchisq qnorm p.adjust
#' @importFrom S4Vectors DataFrame metadata
#' @export
combineVar <- function(..., method=c("z", "fisher", "simes", "berger"), weighted=TRUE) 
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
    output$mean <- .weighted_average_vals(all.means, n.cells, weighted)

    resid.df <- .extract_weightings(all.results, "resid.df")
    for (val in to.average[-1]) { 
        cur.vals <- lapply(all.results, FUN="[[", i=val)
        output[[val]] <- .weighted_average_vals(cur.vals, resid.df, weighted)
    }

    # Combining the one-sided p-values.
    p.combine <- lapply(all.results, FUN="[[", i="p.value")

    method <- match.arg(method)
    if (method=="z") {
        all.z <- lapply(p.combine, FUN=qnorm)
        Z <- .weighted_average_vals(all.z, resid.df, weighted)
        p.final <- pnorm(Z)
    } else if (method=="fisher") {
        all.logp <- lapply(p.combine, FUN=log)
        p.final <- pchisq(-2*Reduce("+", all.logp), df=2*length(all.logp), lower.tail=FALSE)
    } else {
        p.final <- .combine_pvalues(do.call(cbind, p.combine), 
                                    pval.type=ifelse(method=="simes", "any", "all"))
    }

    output <- do.call(DataFrame, output)
    p.final <- unname(p.final)
    output$p.value <- p.final
    output$FDR <- p.adjust(p.final, method="BH")
    rownames(output) <- ref

    # Ensure that you get the same results if you supply only 1 DF as input.
    metadata(output) <- list(num.cells=sum(unlist(n.cells)), 
                             resid.df=sum(unlist(resid.df)))
    return(output)
}

.extract_weightings <- function(tabs, field) {
    output <- vector("list", length(tabs))
    for (x in seq_along(tabs)) { 
        cur.w <- metadata(tabs[[x]])[[field]]
        if (is.null(cur.w)) {
            stop("inputs should come from decomposeVar()")
        }
        output[[x]] <- cur.w
    } 
    names(output) <- names(tabs)
    return(output)
}

.weighted_average_vals <- function(vals, weights, weighted=TRUE) {
    combined <- total.weight <- 0
    for (x in seq_along(vals)) {
        if (weighted) { 
            cur.weights <- weights[[x]]
        } else {
            cur.weights <- 1
        }
        product <- vals[[x]] * cur.weights

        # avoid problems with NA values that have zero weight.
        product[is.na(product) & cur.weights==0] <- 0 

        combined <- combined + product
        total.weight <- total.weight + cur.weights
    }
    combined/total.weight
}
