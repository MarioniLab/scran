#' @importFrom BiocParallel bpmapply SerialParam
#' @importFrom DelayedArray DelayedArray
#' @importFrom DelayedMatrixStats rowVars
#' @importFrom BiocSingular ExactParam
.parallelPCA <- function(x, subset.row=NULL, value=c("pca", "n", "lowrank"), min.rank=5, max.rank=100,
    niters=50, threshold=0.1, BSPARAM=ExactParam(), BPPARAM=SerialParam())
# This performs Horn's parallel analysis to determine the number of PCs
# to retain, by randomizing each row and repeating the PCA to obtain
# an estimate of the mean variance explained per PC under a random model.
#
# written by Aaron Lun
# created 27 March 2018
{
    .Deprecated(old="parallelPCA", new="PCAtools::parallelPCA")

    # Subsetting the matrix.
    x0 <- x
    if (!is.null(subset.row)) {
        subset.row <- .subset_to_index(subset.row, x, byrow=TRUE)
        x <- x[subset.row,,drop=FALSE]
    }
    y <- t(x)
    
    # Running the PCA function once.
    value <- match.arg(value)
    svd.out <- .centered_SVD(y, max.rank, keep.left=(value!="n"), keep.right=(value=="lowrank"), 
        BSPARAM=BSPARAM, BPPARAM=BPPARAM)
    original.d2 <- svd.out$d^2

    # Running it once, and then multiple times after permutation.
    pcg.states <- .setup_pcg_state(niters)
    permuted <- bpmapply(FUN=.parallel_PA, max.rank=rep(max.rank, niters), 
        seed=pcg.states$seeds[[1]], stream=pcg.states$streams[[1]],
        MoreArgs=list(y=y, BSPARAM=BSPARAM), BPPARAM=BPPARAM, SIMPLIFY=FALSE, USE.NAMES=FALSE)
    permutations <- do.call(cbind, permuted)

    # Figuring out where the original drops to "within range" of permuted.
    prop <- rowMeans(permutations >= original.d2)
    above <- prop > threshold
    if (!any(above)) {
        npcs <- length(above) 
    } else {
        npcs <- min(which(above)) - 1L
    }
    npcs <- .keep_rank_in_range(npcs, min.rank, length(original.d2))

    # Collating the return value.
    out.val <- switch(value, 
        n=npcs,
        pca=.svd_to_pca(svd.out, npcs),
        lowrank=.svd_to_lowrank(svd.out, npcs, x0, subset.row)
    )

    var.exp <- original.d2 / (ncol(x) - 1)
    all.var <- sum(rowVars(DelayedArray(x)))
    attr(out.val, "percentVar") <- var.exp/all.var
    attr(out.val, "permuted.percentVar") <- t(permutations)/(ncol(x)-1L)/all.var

    return(out.val)
}

.parallel_PA <- function(y, ..., seed, stream, BSPARAM)
# Function for use in bplapply, defined here to automatically take advantage of the scran namespace when using snowParam. 
# We set keep.left=keep.right=FALSE to avoid computing the left/right eigenvectors, which are unnecessary here.
{
    re.y <- shuffle_matrix(y, seed, stream)
    out <- .centered_SVD(re.y, ..., keep.left=FALSE, keep.right=FALSE, BSPARAM=BSPARAM)
    out$d^2
}

#' @export
setGeneric("parallelPCA", function(x, ...) standardGeneric("parallelPCA"))

#' @export
setMethod("parallelPCA", "ANY", .parallelPCA)

#' @export
#' @importFrom SummarizedExperiment assay "assay<-"
#' @importFrom SingleCellExperiment reducedDim isSpike
setMethod("parallelPCA", "SingleCellExperiment", function(x, ..., subset.row=NULL, value=c("pca", "n", "lowrank"), 
        assay.type="logcounts", get.spikes=FALSE, sce.out=TRUE) {

    subset.row <- .SCE_subset_genes(subset.row=subset.row, x=x, get.spikes=get.spikes)
    out <- .parallelPCA(assay(x, i=assay.type), ..., value=value, subset.row=subset.row)

    value <- match.arg(value) 
    if (!sce.out || value=="n") { 
        return(out)
    }

    if (value=="pca"){ 
        reducedDim(x, "PCA") <- out
    } else if (value=="lowrank") {
        if (!get.spikes) {
            out[isSpike(x),] <- 0
        }
        assay(x, i="lowrank") <- out
    }
    return(x)
})

