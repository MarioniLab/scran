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

    to.average <- c("mean", "total", "bio", "tech")
    output <- list()
    for (val in to.average) { 
        combined <- 0
        for (x in all.results) {
            combined <- combined + x[,val]
        }
        combined <- combined/length(all.results)
        output[[val]] <- combined
    }

    # Combining the p-values.
    p.combine <- list()
    for (i in seq_along(all.results)) {
        p.combine[[i]] <- all.results[[i]]$p.value
    }
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
