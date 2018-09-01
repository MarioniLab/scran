#' @importFrom kmknn findKNN
#' @importFrom methods is
#' @importClassesFrom Matrix dgCMatrix
.quick_sum_factors_per_block <- function(x, k=20, d=50, approximate=FALSE, irlba.args=list(), min.mean=1, subset.row=NULL)
# Implements a much faster method based on local averages to compute size factors.
# Avoids the need for explicit clustering outside of the algorithm.
# 
# written by Aaron Lun
# created 31 August 2018
{
    if (!is.null(subset.row)) {
        x <- x[subset.row,,drop=FALSE]
    }
    if (!is(x, "dgCMatrix")) { # Avoid costly re-reading from file for file-backed matrices.
        x <- as.matrix(x) 
    }
    pcs <- .libnorm_PCA(x, max.rank=d, approximate=approximate, extra.args=irlba.args)
    nn.out <- findKNN(pcs, k=k)
    
    # Choosing the densest cell to be the reference.
    last <- ncol(nn.out$distance)
    if (last==0L) {
        stop("no neighbors available for size factor estimation")
    }
    ref.cell <- which.min(nn.out$distance[,last])

    .quick_sum_cpp_wrapper(x, nn.out$index, nn.out$distance, ref.cell, min.mean=min.mean, ndist=3)
}

#' @importFrom scater librarySizeFactors normalizeCounts
.libnorm_PCA <- function(x, max.rank=50, approximate=FALSE, extra.args=list())
# Performs PCA on the expression profiles after library size normalization.
{
    sf <- librarySizeFactors(x)
    y <- normalizeCounts(x, size_factors=sf, return_log=TRUE)
    SVD <- .centered_SVD(t(y), max.rank=max.rank, approximate=approximate, extra.args=extra.args, keep.right=FALSE)
    .svd_to_pca(SVD, ncomp=max.rank)
}

.quick_sum_cpp_wrapper <- function(x, index, distance, ref.cell, min.mean=NULL, ndist=3) 
# Wrapper to the C++ function for easier calling in the tests.
{
    out <- .Call(scran:::cxx_quick_sum_factors, x, t(index-1L), t(distance), ref.cell - 1L, min.mean, ndist)
    names(out) <- c("sf", "ref")
    out
}

.quick_sum_bpl_wrapper <- function(x, chosen, ...)
# Ensure that namespace is passed along in BiocParallel calls.
{
    .quick_sum_factors_per_block(x[,chosen,drop=FALSE], ...)
}

#' @importFrom BiocParallel SerialParam bpmapply
#' @importFrom stats median
.quickSumFactors <- function(x, k=20, d=50, approximate=FALSE, irlba.args=list(), subset.row=NULL, min.mean=1, block=NULL, BPPARAM=SerialParam()) 
# Parallelizes the size factor calculation across blocks.
{
    all.args <- list(x=x, k=k, d=d, approximate=approximate, irlba.args=irlba.args, subset.row=subset.row, min.mean=min.mean)
    if (is.null(block)) { 
        all.norm <- do.call(.quick_sum_factors_per_block, all.args)
        sf <- all.norm[[1]]
        names(sf) <- colnames(x)
        return(sf/mean(sf))
    }

    indices <- split(seq_along(block), block)
    all.norm <- bpmapply(FUN=.quick_norm_bpl_wrapper, chosen=indices, MoreArgs=all.args, SIMPLIFY=FALSE, USE.NAMES=FALSE, BPPARAM=BPPARAM)

    # Choosing a reference block from the within-block references.
    # This is done by picking the block in the middle of the first PC.
    all.ref <- do.call(rbind, lapply(all.norm, "[[", i="ref"))
    ref.pcs <- .libnorm_PCA(all.ref, max.rank=1, approximate=approximate, extra.args=irlba.args)
    ref.block <- which.min(abs(ref.pcs - median(ref.pcs)))
    true.ref <- all.norm[[ref.block]]$ref
    true.lib <- sum(true.ref)

    # Scaling all size factors to the new reference.
    all.output <- numeric(ncol(x))    
    for (i in seq_along(indices)) {
        cur.sf <- all.norm[[i]]$sf
        cur.ref <- all.norm[[i]]$ref
        cur.lib <- sum(cur.ref)

        ab <- (cur.ref/cur.lib + true.ref/true.lib)/2 * (cur.lib + true.lib)/2
        keep <- ab >= min.mean
        ratios <- cur.ref / true.ref
        all.output[indices[[i]]] <- median(ratios[keep], na.rm=TRUE) * cur.sf
    }

    names(all.output) <- colnames(x)
    return(all.output/mean(all.output))
}

#############################################################
# S4 method definitions.
#############################################################

#' @export
setGeneric("quickSumFactors", function(x, ...) standardGeneric("quickSumFactors"))

#' @export
setMethod("quickSumFactors", "ANY", .quickSumFactors)

#' @importFrom SummarizedExperiment assay 
#' @importFrom BiocGenerics "sizeFactors<-"
#' @export
setMethod("quickSumFactors", "SingleCellExperiment", function(x, ..., min.mean=1, subset.row=NULL, assay.type="counts", get.spikes=FALSE, sf.out=FALSE) 
{ 
    subset.row <- .SCE_subset_genes(subset.row=subset.row, x=x, get.spikes=get.spikes)
    sf <- .quickSumFactors(assay(x, i=assay.type), subset.row=subset.row, min.mean=min.mean, ...) 
    if (sf.out) { 
        return(sf) 
    }
    sizeFactors(x) <- sf
    x
}) 
