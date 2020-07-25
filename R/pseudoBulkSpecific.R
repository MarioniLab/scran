#' Label-specific pseudo-bulk DE
#'
#' Detect label-specific DE genes in a pseudo-bulk analysis,
#' by testing whether the log-fold change is more extreme than the average log-fold change of other labels.
#' 
#' @inheritParams pseudoBulkDGE
#' @param ... For the generic, further arguments to pass to individual methods.
#' 
#' For the ANY method, further arguments to pass to \code{\link{pseudoBulkDGE}}.
#'
#' For the SummarizedExperiment method, further arguments to pass to the ANY method.
#' @param average String specifying the method to compute the average log-fold change of all other labels.
#' @param reference A List containing the (unsorted) output of \code{\link{pseudoBulkDGE}}.
#' This can be supplied to avoid redundant calculations but is automatically computed if \code{NULL}.
#'
#' @details
#' This function implements a quick and dirty method for detecting label-specific DE genes.
#' For a given label and gene, the null hypothesis is that the log-fold change lies between zero
#' and the average log-fold change for that gene across all other labels.
#' This is tested by testing each gene against the two extremes and taking the larger of the two p-values, 
#' as well as setting the p-value to 1 if the log-fold change lies between the extremes.
#'
#' @return
#' A \linkS4class{List} of \linkS4class{DataFrame}s where each DataFrame contains DE statistics for one label.
#' This is equivalent to the output of \code{\link{pseudoBulkDGE}};
#' if \code{reference} is supplied, most of the statistics will be identical to those reported there.
#'
#' The main differences are that the p-values and FDR values are changed.
#' Each DataFrame also has an \code{OtherAverage} field containing the average log-fold change across all other labels.
#'
#' @author Aaron Lun
#'
#' @seealso
#' \code{\link{pseudoBulkDGE}}, for the underlying function that does all the heavy lifting.
#'
#' @examples
#' set.seed(10000)
#' library(scuttle)
#' sce <- mockSCE(ncells=1000)
#' sce$samples <- gl(8, 125) # Pretending we have 8 samples.
#'
#' # Making up some clusters.
#' sce <- logNormCounts(sce)
#' clusters <- kmeans(t(logcounts(sce)), centers=3)$cluster
#'
#' # Creating a set of pseudo-bulk profiles:
#' info <- DataFrame(sample=sce$samples, cluster=clusters)
#' pseudo <- sumCountsAcrossCells(sce, info)
#'
#' # Making up an experimental design for our 8 samples
#' # and adding a common DEG across all labels.
#' pseudo$DRUG <- gl(2,4)[pseudo$sample]
#' assay(pseudo)[1,pseudo$DRUG==1] <- assay(pseudo)[1,pseudo$DRUG==1] * 10 
#' 
#' # Label-specific analysis (note behavior of the first gene).
#' out <- pseudoBulkSpecific(pseudo, 
#'    label=pseudo$cluster,
#'    condition=pseudo$DRUG,
#'    design=~DRUG,
#'    coef="DRUG2"
#' )
#'
#' out[[1]]
#'
#' @export
#' @name pseudoBulkSpecific
NULL

#' @importFrom DelayedMatrixStats rowMedians
#' @importFrom Matrix rowMeans
.pseudo_bulk_specific <- function(x, label, condition, ..., method=c("edgeR", "voom"),
    sorted=FALSE, average=c("median", "mean"), reference=NULL) 
{
    method <- match.arg(method)
    if (is.null(reference)) {
        reference <- pseudoBulkDGE(x, label=label, condition=condition, ..., method=method, sorted=FALSE)
    }

    all.lfc <- lapply(reference, "[[", i="logFC")
    null.lfc.list <- list()
    average <- match.arg(average)

    for (l in names(reference)) {
        other.lfc <- all.lfc[names(all.lfc)!=l]
        other.lfc <- do.call(cbind, other.lfc)

        if (average=="median") {
            null.lfc.list[[l]] <- rowMedians(other.lfc, na.rm=TRUE)
        } else {
            null.lfc.list[[l]] <- rowMeans(other.lfc, na.rm=TRUE)
        }
    }

    alt <- .pseudo_bulk_dge(x=x, label=label, condition=condition,
        null.lfc.list=null.lfc.list, ..., sorted=FALSE, include.intermediates=FALSE)

    for (l in names(reference)) {
        R <- reference[[l]]
        A <- alt[[l]]
        nulls <- null.lfc.list[[l]]

        in.between <- sign(R$logFC) == sign(nulls) & abs(R$logFC) <= abs(nulls)
        p <- ifelse(in.between, 1, pmax(A$PValue, R$PValue))
        fdr <- p.adjust(p, method="BH")

        if (method=="edgeR"){ 
            R$PValue <- p
            R$FDR <- fdr
        } else {
            R$P.Value <- p
            R$adj.P.Val <- fdr
        }

        R$OtherAverage <- nulls

        if (sorted) {
            R <- R[order(p),,drop=FALSE]
        }
        reference[[l]] <- R
    }

    reference
}

#' @export
#' @rdname pseudoBulkSpecific
setGeneric("pseudoBulkSpecific", function(x, ...) standardGeneric("pseudoBulkSpecific"))

#' @export
#' @rdname pseudoBulkSpecific
setMethod("pseudoBulkSpecific", "ANY", .pseudo_bulk_specific)

#' @export
#' @rdname pseudoBulkSpecific
#' @importFrom SummarizedExperiment assay colData
setMethod("pseudoBulkSpecific", "SummarizedExperiment", function(x, ..., assay.type=1) {
    .pseudo_bulk_specific(assay(x, assay.type), col.data=colData(x), ...)
})
