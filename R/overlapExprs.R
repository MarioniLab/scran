#' @importFrom BiocParallel SerialParam bplapply
.overlapExprs <- function(x, groups, gene.names=rownames(x), block=NULL, pval.type=c("any", "all"),
    direction=c("any", "up", "down"), tol=1e-8, log.p=FALSE, full.stats=FALSE, subset.row=NULL, BPPARAM=SerialParam()) 
# Computes the gene-specific overlap in expression profiles between two groups of cells.
# This aims to determine whether two distributions of expression values are well-separated.    
# 
# written by Aaron Lun
# created 17 April 2017
{
    fit <- pairwiseWilcox(x, groups, block=block, direction=direction, gene.names=gene.names, log.p=TRUE, subset.row=subset.row, BPPARAM=BPPARAM)
    combineMarkers(fit$statistics, fit$pairs, pval.type=pval.type, log.p.in=TRUE, log.p.out=log.p, full.stats=full.stats, pval.field="log.p.value", effect.field="AUC")
}

###########################################################
# S4 method definitions
###########################################################

#' @export
setGeneric("overlapExprs", function(x, ...) standardGeneric("overlapExprs"))

#' @export
setMethod("overlapExprs", "ANY", .overlapExprs)

#' @importFrom SummarizedExperiment assay
#' @export
setMethod("overlapExprs", "SingleCellExperiment", function(x, ..., subset.row=NULL, assay.type="logcounts", get.spikes=FALSE) {
    subset.row <- .SCE_subset_genes(subset.row=subset.row, x=x, get.spikes=get.spikes)
    .overlapExprs(assay(x, i=assay.type), ..., subset.row=subset.row)
})                                 

