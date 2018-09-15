.findMarkers <- function(x, clusters, gene.names=rownames(x), block=NULL, design=NULL, 
    pval.type=c("any", "all"), direction=c("any", "up", "down"), 
    lfc=0, log.p=FALSE, full.stats=FALSE, subset.row=NULL)
# Uses limma to find the markers that are differentially expressed between clusters,
# given a log-expression matrix and some blocking factors in 'design' or 'block'.
#
# written by Aaron Lun
# created 22 March 2017
{
    fit <- pairwiseTTests(x, clusters, block=block, design=design, direction=direction, lfc=lfc, gene.names=gene.names, log.p=TRUE, subset.row=subset.row)
    combineMarkers(fit$statistics, fit$pairs, pval.type=pval.type, log.p.in=TRUE, log.p.out=log.p, full.stats=full.stats, pval.field="log.p.value")
}

#' @export
setGeneric("findMarkers", function(x, ...) standardGeneric("findMarkers"))

#' @export
setMethod("findMarkers", "ANY", .findMarkers)

#' @importFrom SummarizedExperiment assay
#' @export
setMethod("findMarkers", "SingleCellExperiment", function(x, ..., subset.row=NULL, assay.type="logcounts", get.spikes=FALSE) {
    subset.row <- .SCE_subset_genes(subset.row=subset.row, x=x, get.spikes=get.spikes)
    .findMarkers(assay(x, i=assay.type), ..., subset.row=subset.row)
})
