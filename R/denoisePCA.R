.denoisePCA <- function(x, technical, design=NULL, subset.row=NULL,
                        value=c("pca", "n", "lowrank"), min.rank=5, max.rank=100)
# Performs PCA and chooses the number of PCs to keep based on the technical noise.
# This is done on the residuals if a design matrix is supplied.
#
# written by Aaron Lun
# created 13 March 2017    
# last modified 15 July 2017
{
    subset.row <- .subset_to_index(subset.row, x, byrow=TRUE)
    x <- x[subset.row,] # Might as well, need to do PCA on the subsetted matrix anyway.
    subset.row <- seq_len(nrow(x))
    all.means <- rowMeans(x)

    if (!is.null(design)) { 
        checked <- .make_var_defaults(x, fit=NULL, design=design)
        design <- checked$design
        QR <- .ranksafe_qr(design)
        
        # Computing residuals; don't set a lower bound.
        rx <- .calc_residuals_wt_zeroes(x, QR=QR, subset.row=subset.row, lower.bound=NA) 

        # Rescaling residuals so that the variance is unbiased.
        # This is necessary because variance of residuals is underestimated.
        xout <- .Call(cxx_estimate_variance, QR$qr, QR$qraux, x, subset.row - 1L)
        xvar <- xout[[2]]
        rvar <- apply(rx, 1, var)
       
        # Replacing 'x' with the scaled residuals (these shoud have a mean of zero,
        # see http://math.stackexchange.com/questions/494181/ for a good explanation).
        x <- rx * sqrt(xvar/rvar)
    }

    # Computing the technical variance sum.
    if (is.function(technical)) { 
        technical <- sum(technical(all.means))
    } else if (is.numeric(technical)) { 
        if (is.null(rownames(x))) { 
            stop("rows of 'x' should be named with gene names")
        }
        technical <- sum(technical[rownames(x)])
        if (is.na(technical)) {
            stop("missing gene names in 'technical'")
        }
    } else {
        stop("'technical' should be a function or a scalar")
    }

    # Performing SVD to get the variance of each PC, and choosing the number of PCs to keep.
    centers <- rowMeans(x)
    y <- t(x - centers)
    svd.out <- svd(y, nu=0, nv=0)
    var.exp <- svd.out$d^2/(ncol(x) - 1)

    to.keep <- .get_npcs_to_keep(var.exp, technical)
    to.keep <- min(max(to.keep, min.rank), max.rank)
    
    # Figuring out what value to return; the number of PCs, the PCs themselves, or a denoised low-rank matrix.
    value <- match.arg(value)
    if (value=="n") {
        return(to.keep)
    } else if (value=="pca") {
        pc.out <- prcomp(y, rank.=to.keep, scale.=FALSE, center=FALSE)
        return(pc.out$x)
    } else if (value=="lowrank") {
        more.svd <- La.svd(y, nu=to.keep, nv=to.keep)
        denoised <- more.svd$u %*% (more.svd$d[seq_len(to.keep)] * more.svd$vt) 
        denoised <- t(denoised) + centers
        dimnames(denoised) <- dimnames(x)
        return(denoised)
    }
} 

.get_npcs_to_keep <- function(var.exp, tech.var) 
# Discarding PCs until we get rid of as much technical noise as possible
# while preserving the biological signal. This is done by assuming that 
# the biological signal is fully contained in earlier PCs, such that we 
# discard the later PCs until we account for 'tech.var'.
{
    npcs <- length(var.exp)
    flipped.var.exp <- rev(var.exp)
    estimated.contrib <- cumsum(flipped.var.exp) + flipped.var.exp * (npcs:1 - 1L)
    estimated.contrib <- rev(estimated.contrib)

    below.noise <- tech.var > estimated.contrib
    if (any(below.noise)) { 
        to.keep <- min(which(below.noise))
    } else {
        to.keep <- npcs
    }
    return(to.keep)
}

setGeneric("denoisePCA", function(x, ...) standardGeneric("denoisePCA"))

setMethod("denoisePCA", "matrix", .denoisePCA)

setMethod("denoisePCA", "SCESet", function(x, ..., subset.row=NULL, value=c("pca", "n", "lowrank"), 
                                           assay="exprs", get.spikes=FALSE) {
    if (is.null(subset.row)) {
        subset.row <- .spike_subset(x, get.spikes)
    }
    out <- .denoisePCA(assayDataElement(x, assay), ..., value=value, subset.row=subset.row)

    value <- match.arg(value) 
    if (value=="pca"){ 
        reducedDimension(x) <- out
    } else if (value=="n") {
        ; # will put this into metadata in the future.
    } else if (value=="lowrank") {
        original <- assayDataElement(x, assay)
        subset.row <- .subset_to_index(subset.row, original, byrow=TRUE)
        output <- matrix(NA_real_, nrow(original), ncol(original))
        output[subset.row,] <- out
        assayDataElement(x, "lowrank") <- output
    }
    return(x)
})

