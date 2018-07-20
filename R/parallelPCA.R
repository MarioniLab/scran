#' @importFrom BiocParallel bplapply SerialParam
#' @importFrom DelayedArray DelayedArray
#' @importFrom DelayedMatrixStats rowVars
.parallelPCA <- function(x, subset.row=NULL, scale=NULL, value=c("pca", "n", "lowrank"), min.rank=5, max.rank=100,
                         niters=50, threshold=0.1, approximate=FALSE, irlba.args=list(), BPPARAM=SerialParam())
# This performs Horn's parallel analysis to determine the number of PCs
# to retain, by randomizing each row and repeating the PCA to obtain
# an estimate of the mean variance explained per PC under a random model.
#
# written by Aaron Lun
# created 27 March 2018
{
    x0 <- x
    scale0 <- scale

    # Subsetting and scaling the matrix.
    if (!is.null(subset.row)) {
        subset.row <- .subset_to_index(subset.row, x, byrow=TRUE)
        x <- x[subset.row,,drop=FALSE]
        scale <- scale[subset.row]
    }

    if (!is.null(scale)) {
        x <- x * scale
    }
    
    # Setting up the PCA function and its arguments.
    value <- match.arg(value)
    args <- list(y=t(x), max.rank=max.rank, value=value)
    if (approximate) {
        svdfun <- .irlba_svd
        args <- c(args, irlba.args)
    } else {
        svdfun <- .full_svd
    }

    # Running it once, and then multiple times after permutation.
    original <- do.call(svdfun, args)
    original.d2 <- original$d^2
    permuted <- bplapply(seq_len(niters), FUN=.parallel_PA, svdfun=svdfun, args=args, BPPARAM=BPPARAM)
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
    out.val <- .convert_to_output(original, npcs, value, x0, scale0, subset.row)

    var.exp <- original.d2 / (ncol(x) - 1)
    all.var <- sum(rowVars(DelayedArray(x)))
    attr(out.val, "percentVar") <- var.exp/all.var
    attr(out.val, "permuted.percentVar") <- t(permutations)/(ncol(x)-1L)/all.var

    return(out.val)
}

.parallel_PA <- function(svdfun, args, ...) 
# Function for use in bplapply, defined here to automatically take 
# advantage of the scran namespace when using snowParam. We set
# value='n' as we don't really want anything but the singular values here.
{
    args$value <- "n"
    args$y <- .Call(cxx_shuffle_matrix, args$y)
    out <- do.call(svdfun, args)
    return(out$d^2)
}

#' @export
setGeneric("parallelPCA", function(x, ...) standardGeneric("parallelPCA"))

#' @export
setMethod("parallelPCA", "ANY", .parallelPCA)

#' @importFrom SummarizedExperiment assay "assay<-"
#' @importFrom SingleCellExperiment reducedDim isSpike
#' @export
setMethod("parallelPCA", "SingleCellExperiment", 
          function(x, ..., subset.row=NULL, value=c("pca", "n", "lowrank"), 
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

