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
#' @param missing.as.zero Logical scalar indicating whether missing log-fold changes should be set to zero.
#' @param reference A List containing the (unsorted) output of \code{\link{pseudoBulkDGE}}.
#' This can be supplied to avoid redundant calculations but is automatically computed if \code{NULL}.
#'
#' @details
#' This function implements a quick and dirty method for detecting label-specific DE genes.
#' For a given label and gene, the null hypothesis is that the log-fold change lies between zero
#' and the average log-fold change for that gene across all other labels.
#' Genes that reject this null either have log-fold changes in the opposite sign 
#' or are significantly more extreme than the average.
#'
#' To implement this, we test each gene against the two extremes and taking the larger of the two p-values.
#' The p-value is set to 1 if the log-fold change lies between the extremes.
#' This is somewhat similar to how \code{\link{treat}} might behave if the null interval was not centered at zero;
#' however, our approach is more conservative than the \code{\link{treat}} as the p-value calculations are not quite correct.
#'
#' It is worth stressing that there are no guarantees that the DE genes detected in this manner are truly label-specific.
#' For any label and DEG, there may be one or more labels with stronger log-fold changes,
#' but the average may be pulled towards zero by other labels with weaker or opposing effects.
#' The use of the average is analogous to recommendations in the \pkg{edgeR} user's guide for testing against multiple groups.
#' However, a more stringent selection can be achieved by applying gates on \code{\link{decideTestsPerLabel}}.
#'
#' Note that, if \code{lfc} is specified in the arguments to \code{\link{pseudoBulkDGE}},
#' the null interval is expanded in both directions by the specified value.
#' 
#' @section Computing the average:
#' The average log-fold change for each gene is computed by taking the median or mean (depending on \code{average}) 
#' of the corresponding log-fold changes in each of the DE analyses for the other labels.
#' We use the median by default as this means that at least half of all other labels should have weaker or opposing effects.
#'
#' By default, low-abundance genes that were filtered out in a comparison do not contribute to the average.
#' Any log-fold changes that could be computed from them are considered to be too unstable.
#' If the gene is filtered out in all other labels, the average is set to zero for testing but is reported as \code{NA}. 
#'
#' If \code{missing.as.zero=TRUE}, the log-fold changes for all filtered genes are set to zero.
#' This is useful when a gene is only expressed in the subset of labels and is consistently DEG in each comparison of the subset.
#' Testing against the average computed from only those labels in the subset would fail to detect this DEG as subset-specific.
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
.pseudo_bulk_specific <- function(x, label, condition=NULL, ..., method=c("edgeR", "voom"),
    sorted=FALSE, average=c("median", "mean"), missing.as.zero=FALSE, reference=NULL) 
{
    method <- match.arg(method)
    if (is.null(reference)) {
        reference <- pseudoBulkDGE(x, label=label, condition=condition, ..., method=method, sorted=FALSE)
    }

    all.lfc <- lapply(reference, "[[", i="logFC")
    if (any(vapply(all.lfc, is.null, FALSE))) {
        stop("ANOVA-like contrasts cannot be specified with non-zero 'null.lfc'")
    }

    null.lfc.list <- lost <- list()
    average <- match.arg(average)

    for (l in names(reference)) {
        other.lfc <- all.lfc[names(all.lfc)!=l]
        other.lfc <- do.call(cbind, other.lfc)

        if (missing.as.zero) {
            other.lfc[is.na(other.lfc)] <- 0
        }

        if (average=="median") {
            ave.other <- rowMedians(other.lfc, na.rm=TRUE)
        } else {
            ave.other <- rowMeans(other.lfc, na.rm=TRUE)
        }

        zeroed <- is.na(ave.other)
        ave.other[zeroed] <- 0
        lost[[l]] <- zeroed
        null.lfc.list[[l]] <- ave.other
    }

    alt <- .pseudo_bulk_dge(x=x, label=label, condition=condition, method=method,
        null.lfc.list=null.lfc.list, ..., sorted=FALSE, include.intermediates=FALSE)

    if (method=="edgeR"){
        pname <- "PValue"
        fname <- "FDR"
    } else {
        pname <- "P.Value"
        fname <- "adj.P.Val"
    }

    for (l in names(reference)) {
        R <- reference[[l]]
        A <- alt[[l]]
        nulls <- null.lfc.list[[l]]

        in.between <- sign(R$logFC) == sign(nulls) & abs(R$logFC) <= abs(nulls)
        p <- ifelse(in.between, 1, pmax(A[[pname]], R[[pname]]))
        R[[pname]] <- p
        R[[fname]] <- p.adjust(p, method="BH")

        R$OtherAverage <- ifelse(lost[[l]], NA_real_, nulls)

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
