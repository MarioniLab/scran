.denoisePCA <- function(x, technical, design=NULL, subset.row=NULL,
                        value=c("pca", "n", "lowrank"), min.rank=5, max.rank=100, 
                        approximate=FALSE, rand.seed=1000, irlba.args=list())
# Performs PCA and chooses the number of PCs to keep based on the technical noise.
# This is done on the residuals if a design matrix is supplied.
#
# written by Aaron Lun
# created 13 March 2017    
# last modified 15 November 2017
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
        centering <- numeric(nrow(y))
    } else {
        y <- x[use.rows,,drop=FALSE] 
        centering <- all.means # all.means is already subsetted, remember.
    }

    # Checking various other arguments.
    value <- match.arg(value)
    min.rank <- max(1L, min.rank)

    # Switching to IRLBA if an approximation is requested.
    if (approximate) {
        if (!is.na(rand.seed)) { 
            set.seed(rand.seed)
        }
        out.val <- .approximate_denoisePCA(y=y, min.rank=min.rank, max.rank=max.rank,
            value=value, technical=technical, centering=centering, 
            all.var=all.var, irlba.args=irlba.args)
    } else {
        out.val <- .exact_denoisePCA(y=y, min.rank=min.rank, max.rank=max.rank, 
            value=value, technical=technical, centering=centering)
    }
    var.exp <- out.val$variance
    out.val <- out.val$output
    
    # Computing a low-rank approximation here.
    if (value=="lowrank") { 
        denoised <- out.val$U %*% (out.val$D * out.val$Vt) 
        denoised <- t(denoised) + all.means
        
        output <- matrix(0, nrow=nrow(x), ncol=ncol(x))
        dimnames(output) <- dimnames(x)
        output[use.rows,] <- denoised

        # The idea is that after our SVD, we have X=UDV' where each column of X is a gene. 
        # Leftover genes are new columns in X, which are projected on the space of U by doing U'X.
        # This can be treated as new columns in DV', which can be multiplied by U to give denoised values.
        # I've done a lot of implicit transpositions here, hence the code does not tightly follow the logic above.
        leftovers <- seq_len(nrow(x))[-use.rows]
        if (length(leftovers)) { 
            left.x <- x[leftovers,,drop=FALSE]
            left.means <- calcAverage(left.x, size.factors=rep(1, ncol(x)))
            left.x <- left.x - left.means
            new.vals <- tcrossprod(left.x %*% out.val$U, out.val$U)
            new.vals <- new.vals + left.means
            output[leftovers,] <- new.vals
        }

        out.val <- output
    }

    # Adding the percentage of variance explained.
    attr(out.val, "percentVar") <- var.exp/sum(all.var)
    return(out.val)
} 

##################################
# Internal SVD and PCA functions #
##################################

.approximate_denoisePCA <- function(y, min.rank, max.rank, value, technical, centering, all.var, irlba.args) { 
    ncells <- ncol(y)
    max.rank <- min(dim(y)-1L, max.rank)
    nu <- ifelse(value!="n", max.rank, 0L)

    # Allowing more iterations if max.rank is high.
    arg.max <- pmatch(names(irlba.args), "maxit")
    if (all(is.na(arg.max))) { 
        irlba.args$maxit <- max(100, max.rank*10)
    }

    out <- do.call(irlba::irlba, c(list(A=t(y), nu=nu, nv=max.rank, center=centering), irlba.args)) 
    var.exp <- out$d^2/(ncells - 1)
    
    # Assuming all non-computed components were technical, and discarding them for further consideration.
    leftovers <- sum(all.var) - sum(var.exp)
    technical <- max(0, technical - leftovers)
    to.keep <- .get_npcs_to_keep(var.exp, technical)
    to.keep <- min(max(to.keep, min.rank), max.rank)
    
    # Figuring out what value to return; the number of PCs, the PCs themselves, or a denoised low-rank matrix.
    if (value=="n") {
        out.val <- to.keep
    } else if (value=="pca") {
        ix <- seq_len(to.keep)
        out.val <- sweep(out$u[,ix,drop=FALSE], 2, out$d[ix], FUN="*")
    } else if (value=="lowrank") {
        ix <- seq_len(to.keep)
        out.val <- list(U=out$u[,ix,drop=FALSE], D=out$d[ix], Vt=t(out$v[,ix,drop=FALSE]))
    }
    return(list(output=out.val, variance=var.exp))
}

.exact_denoisePCA <- function(y, min.rank, max.rank, value, technical, centering) { 
    ncells <- ncol(y)
    max.rank <- min(dim(y), max.rank)

    # Centering the matrix and coercing it to a dense representation.
    y <- t(y - centering)
    y <- as.matrix(y)

    # Performing SVD to get the variance of each PC, and choosing the number of PCs to keep.
    svd.out <- svd(y, nu=0, nv=0)
    var.exp <- svd.out$d^2/(ncells - 1)
    to.keep <- .get_npcs_to_keep(var.exp, technical)
    to.keep <- min(max(to.keep, min.rank), max.rank)
    
    # Figuring out what value to return; the number of PCs, the PCs themselves, or a denoised low-rank matrix.
    if (value=="n") {
        out.val <- to.keep
    } else if (value=="pca") {
        pc.out <- prcomp(y, rank.=to.keep, scale.=FALSE, center=FALSE)
        out.val <- pc.out$x
    } else if (value=="lowrank") {
        more.svd <- La.svd(y, nu=to.keep, nv=to.keep)
        out.val <- list(U=more.svd$u, D=more.svd$d[seq_len(to.keep)], Vt=more.svd$vt)
    }
    return(list(output=out.val, variance=var.exp))
}

.get_npcs_to_keep <- function(var.exp, tech.var) 
# Discarding PCs until we get rid of as much technical noise as possible
# while preserving the biological signal. This is done by assuming that 
# the biological signal is fully contained in earlier PCs, such that we 
# discard the later PCs until we account for 'tech.var'.
{
    npcs <- length(var.exp)
    flipped.var.exp <- rev(var.exp)
    estimated.contrib <- cumsum(flipped.var.exp) 

    above.noise <- estimated.contrib > tech.var
    if (any(above.noise)) { 
        to.keep <- npcs - min(which(above.noise)) + 1L
    } else {
        to.keep <- npcs
    }
    return(to.keep)
}

##############################
# S4 method definitions here #
##############################

setGeneric("denoisePCA", function(x, ...) standardGeneric("denoisePCA"))

setMethod("denoisePCA", "ANY", .denoisePCA)

setMethod("denoisePCA", "SingleCellExperiment", 
          function(x, ..., subset.row=NULL, value=c("pca", "n", "lowrank"), 
                   assay.type="logcounts", get.spikes=FALSE, sce.out=TRUE) {

    subset.row <- .SCE_subset_genes(subset.row=subset.row, x=x, get.spikes=get.spikes)
    out <- .denoisePCA(assay(x, i=assay.type), ..., value=value, subset.row=subset.row)

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

