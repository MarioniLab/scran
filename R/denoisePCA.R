.denoisePCA <- function(x, technical, design=NULL, subset.row=NULL,
                        value=c("pca", "n", "lowrank"), min.rank=5, max.rank=100, 
                        preserve.dim=FALSE, approximate=FALSE, rand.seed=1000)
# Performs PCA and chooses the number of PCs to keep based on the technical noise.
# This is done on the residuals if a design matrix is supplied.
#
# written by Aaron Lun
# created 13 March 2017    
# last modified 17 August 2017
{
    subset.row <- .subset_to_index(subset.row, x, byrow=TRUE)
    checked <- .make_var_defaults(x, fit=NULL, design=design)
    QR <- .ranksafe_qr(checked$design)
    stats <- .Call(cxx_fit_linear_model, QR$qr, QR$qraux, x, subset.row - 1L, FALSE)
    all.means <- stats[[1]]
    all.var <- stats[[2]]

    # Filtering out genes with negative biological components.
    tech.var <- technical(all.means)
    keep <- all.var > tech.var
    use.rows <- subset.row[keep]
    all.means <- all.means[keep]
    all.var <- all.var[keep]
    tech.var <- tech.var[keep]
    technical <- sum(tech.var)

    if (!is.null(design)) {  
        # Computing residuals; don't set a lower bound.
        # Note that this function implicitly subsets by subset.row.
        rx <- .calc_residuals_wt_zeroes(x, QR=QR, subset.row=use.rows, lower.bound=NA) 

        # Rescaling residuals so that the variance is unbiased.
        # This is necessary because variance of residuals is underestimated.
        rvar <- apply(rx, 1, var)
       
        # Replacing 'x' with the scaled residuals (these should already have a mean of zero,
        # see http://math.stackexchange.com/questions/494181/ for a good explanation).
        y <- rx * sqrt(all.var/rvar)
    } else {
        y <- x[use.rows,,drop=FALSE] - all.means 
    }

    # Checking various other arguments.
    value <- match.arg(value)
    min.rank <- max(1L, min.rank)
    ncells <- ncol(x)
    max.rank <- min(ncells, max.rank)
    y <- t(y)

    # Switching to IRLBA if an approximation is requested.
    if (approximate){
        if (!is.na(rand.seed)) { 
            set.seed(rand.seed)
        }
        max.rank <- min(dim(y)-1L, max.rank)
        nu <- ifelse(value!="n", max.rank, 0L)
        out <- irlba::irlba(y, nu=nu, nv=max.rank, # center= seems to be broken.
                            maxit=max(100, max.rank*10)) # allowing more iterations if max.rank is high.
        var.exp <- out$d^2/(ncells - 1)
        
        # Assuming all non-computed components were technical, and discarding them for further consideration.
        leftovers <- sum(all.var) - sum(var.exp)
        technical <- max(0, technical - leftovers)
        to.keep <- .get_npcs_to_keep(var.exp, technical)
        to.keep <- min(max(to.keep, min.rank), max.rank)
        
        # Figuring out what value to return; the number of PCs, the PCs themselves, or a denoised low-rank matrix.
        if (value=="n") {
            return(to.keep)
        } else if (value=="pca") {
            ix <- seq_len(to.keep)
            return(sweep(out$u[,ix,drop=FALSE], 2, out$d[ix], FUN="*"))
        } else if (value=="lowrank") {
            ix <- seq_len(to.keep)
            denoised <- out$u[,ix,drop=FALSE] %*% (out$d[ix] * t(out$v[,ix,drop=FALSE])) 
            denoised <- t(denoised) + all.means
            return(.restore_dimensions(x, denoised, use.rows, subset.row, preserve.dim=preserve.dim))
        }
    }

    # Performing SVD to get the variance of each PC, and choosing the number of PCs to keep.
    svd.out <- svd(y, nu=0, nv=0)
    var.exp <- svd.out$d^2/(ncells - 1)
    to.keep <- .get_npcs_to_keep(var.exp, technical)
    to.keep <- min(max(to.keep, min.rank), max.rank)
    
    # Figuring out what value to return; the number of PCs, the PCs themselves, or a denoised low-rank matrix.
    if (value=="n") {
        return(to.keep)
    } else if (value=="pca") {
        pc.out <- prcomp(y, rank.=to.keep, scale.=FALSE, center=FALSE)
        return(pc.out$x)
    } else if (value=="lowrank") {
        more.svd <- La.svd(y, nu=to.keep, nv=to.keep)
        denoised <- more.svd$u %*% (more.svd$d[seq_len(to.keep)] * more.svd$vt) 
        denoised <- t(denoised) + all.means
        return(.restore_dimensions(x, denoised, use.rows, subset.row, preserve.dim=preserve.dim))
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

    below.noise <- tech.var > estimated.contrib
    if (any(below.noise)) { 
        to.keep <- npcs - max(which(below.noise))
    } else {
        to.keep <- npcs
    }
    return(to.keep)
}

.restore_dimensions <- function(original, denoised, use.rows, subset.row, preserve.dim=FALSE)
# Returning as a the full matrix (or that subsetted with subset.row)
# where the discarded genes are set to zero.
{
    output <- matrix(0, nrow=nrow(original), ncol=ncol(original))
    dimnames(output) <- dimnames(original)
    if (preserve.dim) { 
        output[use.rows,] <- denoised
    } else {
        output <- output[subset.row,,drop=FALSE]
        output[match(use.rows, subset.row),] <- denoised
    }
    return(output)
}

setGeneric("denoisePCA", function(x, ...) standardGeneric("denoisePCA"))

setMethod("denoisePCA", "ANY", .denoisePCA)

setMethod("denoisePCA", "SingleCellExperiment", 
          function(x, ..., subset.row=NULL, value=c("pca", "n", "lowrank"), 
                   assay.type="logcounts", get.spikes=FALSE, sce.out=TRUE) {

    subset.row <- .SCE_subset_genes(subset.row=subset.row, x=x, get.spikes=get.spikes)
    out <- .denoisePCA(assay(x, i=assay.type), ..., value=value, subset.row=subset.row, preserve.dim=TRUE)

    value <- match.arg(value) 
    if (!sce.out || value=="n") { 
        return(out)
    }

    if (value=="pca"){ 
        reducedDim(x, "PCA") <- out
    } else if (value=="lowrank") {
        assay(x, i="lowrank") <- out
    }
    return(x)
})

