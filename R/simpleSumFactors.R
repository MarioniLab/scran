#' @importFrom BiocNeighbors findKNN
#' @importFrom methods is
#' @importClassesFrom Matrix dgCMatrix
#' @importFrom scater librarySizeFactors normalizeCounts
#' @importFrom BiocGenerics t
.simple_sum_factors_per_block <- function(x, indices=NULL, k=20, trend.args=list(), approximate=FALSE, irlba.args=list(), min.mean=1, subset.row=NULL, BNPARAM=NULL)
# Implements a much faster method based on local averages to compute size factors.
# Avoids the need for explicit clustering outside of the algorithm.
# 
# written by Aaron Lun
# created 31 August 2018
{
    if (!is.null(indices)) {
        x <- x[,indices,drop=FALSE]
    }
    if (!is.null(subset.row)) {
        x <- x[subset.row,,drop=FALSE]
    }
    if (!is(x, "dgCMatrix")) { # Avoid costly re-reading from file for file-backed matrices.
        x <- as.matrix(x) 
    }
    if (ncol(x)==0L) {
        stop("no neighbors available for size factor estimation")
    }

    # A mini-analysis based on library size normalization.
    sf <- librarySizeFactors(x)
    y <- normalizeCounts(x, size_factors=sf, return_log=TRUE)
    fit <- do.call(trendVar, c(list(x=y), trend.args))
    pcs <- denoisePCA(y, technical=fit$trend, approximate=approximate, irlba.args=irlba.args)
    nn.out <- findKNN(pcs, k=k, BNPARAM=BNPARAM)
    
    # Choosing the densest cell to be the reference.
    last <- ncol(nn.out$distance)
    ref.cell <- which.min(nn.out$distance[,last])

    .simple_sum_cpp_wrapper(x, nn.out$index, nn.out$distance, ref.cell, min.mean=min.mean, ndist=3)
}

.simple_sum_cpp_wrapper <- function(x, index, distance, ref.cell, min.mean=NULL, ndist=3) 
# Wrapper to the C++ function for easier calling in the tests.
{
    out <- .Call(cxx_simple_sum_factors, x, t(index-1L), t(distance), ref.cell - 1L, min.mean, ndist)
    names(out) <- c("sf", "ref")
    out
}

#' @importFrom BiocParallel SerialParam bpmapply
#' @importFrom stats median
#' @importFrom scater librarySizeFactors normalizeCounts
#' @importFrom BiocGenerics t
.simpleSumFactors <- function(x, k=20, trend.args=list(), approximate=FALSE, irlba.args=list(), subset.row=NULL, min.mean=1, block=NULL, BNPARAM=NULL, BPPARAM=SerialParam()) 
# Parallelizes the size factor calculation across blocks.
{
    all.args <- list(x=x, k=k, trend.args=trend.args, approximate=approximate, irlba.args=irlba.args, subset.row=subset.row, min.mean=min.mean, BNPARAM=BNPARAM)
    if (is.null(block)) { 
        all.norm <- do.call(.simple_sum_factors_per_block, all.args)
        sf <- all.norm[[1]]
        names(sf) <- colnames(x)
        return(sf/mean(sf))
    }

    indices <- split(seq_along(block), block)
    all.norm <- bpmapply(FUN=.simple_sum_factors_per_block, indices=indices, MoreArgs=all.args, SIMPLIFY=FALSE, USE.NAMES=FALSE, BPPARAM=BPPARAM)

    # Choosing a reference block from the within-block references.
    # This is done by picking the block in the middle of the first PC.
    all.ref <- lapply(all.norm, "[[", i="ref")
    R <- do.call(cbind, all.ref)
    R <- normalizeCounts(R, size_factors=librarySizeFactors(R), return_log=TRUE)
    ref.pcs <- .centered_SVD(t(R), max.rank=1, approximate=approximate, extra.args=irlba.args)$u
    ref.block <- which.min(abs(ref.pcs - median(ref.pcs)))

    # Scaling all size factors to the new reference.
    rescaling.factors <- .rescale_clusters(all.ref, ref.clust=ref.block, min.mean=min.mean, clust.names=names(indices))
    all.output <- numeric(ncol(x))
    for (block in seq_along(all.ref)) {
        all.output[indices[[block]]] <- all.norm[[block]]$sf * rescaling.factors[[block]]
    }

    names(all.output) <- colnames(x)
    return(all.output/mean(all.output))
}

#############################################################
# S4 method definitions.
#############################################################

#' @export
setGeneric("simpleSumFactors", function(x, ...) standardGeneric("simpleSumFactors"))

#' @export
setMethod("simpleSumFactors", "ANY", .simpleSumFactors)

#' @importFrom SummarizedExperiment assay 
#' @importFrom BiocGenerics "sizeFactors<-"
#' @export
setMethod("simpleSumFactors", "SingleCellExperiment", function(x, ..., subset.row=NULL, assay.type="counts", get.spikes=FALSE, sf.out=FALSE) 
{ 
    subset.row <- .SCE_subset_genes(subset.row=subset.row, x=x, get.spikes=get.spikes)
    sf <- .simpleSumFactors(assay(x, i=assay.type), subset.row=subset.row, ...) 
    if (sf.out) { 
        return(sf) 
    }
    sizeFactors(x) <- sf
    x
}) 
