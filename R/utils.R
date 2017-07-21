.spike_subset <- function(x, get.spikes) {
    if (!get.spikes) {
        nokeep <- isSpike(x, warning=FALSE)
        if (!is.null(nokeep) && any(nokeep)) {
            return(!nokeep)
        }
    } 
    return(NULL)
}

.subset_to_index <- function(subset, x, byrow=TRUE) {
    if (byrow) {
        dummy <- seq_len(nrow(x))
        names(dummy) <- rownames(x)
    } else {
        dummy <- seq_len(ncol(x))
        names(dummy) <- colnames(x) 
    }

    if (!is.null(subset)) { 
        dummy <- dummy[subset]
    }
    out <- unname(dummy)
    if (any(is.na(out))) {
        stop("'subset' indices out of range of 'x'")
    }
    return(out)
}

#################################################

.make_var_defaults <- function(x, fit, design) 
# Makes defaults for the trendVar and decomposeVar functions.
{
    if (is.null(design)) { 
       design <- as.matrix(rep(1, ncol(x))) 
    } else if (length(design)==1L && is.na(design)) { 
        design <- fit$design 
    }
    return(list(design=design))
}

.check_centered_SF <- function(x, assay) 
# Checks if 'exprs' was requested, and if it could have been computed from counts,
# If so, then it checks if the size factors are centered.
{
    if (assay=="exprs" && 
        !is.null(suppressWarnings(get_exprs(x, "counts", warning=TRUE))) && 
        !areSizeFactorsCentred(x)) {
        warning("size factors not centred, run 'normalize()' first")
    }
    return(NULL)
}

.prepare_cv2_data <- function(x, spike.type) 
# Prepares data for calculation of CV2.
# In particular, extracting spike-ins and their size factors.    
{
    sf.cell <- sizeFactors(x)
    if (is.null(spike.type) || !is.na(spike.type)) { 
        is.spike <- isSpike(x, type=spike.type)
        if (is.null(spike.type)) { 
            # Get all spikes.
            spike.type <- whichSpike(x)            
        }
        if (!length(spike.type)) { 
            stop("no spike-in sets specified from 'x'")
        }

        # Collecting the size factors for the requested spike-in sets.
        # Check that all spike-in factors are either NULL or identical.
        collected <- NULL
        for (st in seq_along(spike.type)) {
            cur.sf <- suppressWarnings(sizeFactors(x, type=spike.type[st]))
            if (st==1L) {
                collected <- cur.sf
            } else if (!isTRUE(all.equal(collected, cur.sf))) {
                stop("size factors differ between spike-in sets")
            }
        }

        # Otherwise, diverting to the cell-based size factors if all spike-in factors are NULL.
        if (!is.null(collected)) {
            sf.spike <- collected
        } else {
            warning("no spike-in size factors set, using cell-based factors")
            sf.spike <- sf.cell
        }

    } else {
        sf.spike <- sf.cell
        is.spike <- NA
    }

    return(list(is.spike=is.spike, sf.cell=sf.cell, sf.spike=sf.spike))
} 

#################################################

.worker_assign <- function(njobs, BPPARAM)
# Assigns jobs to workers.
{
    ncores <- bpworkers(BPPARAM)
    starting <- as.integer(seq(from=1, to=njobs+1, length.out=ncores+1))
    jobsize <- diff(starting)
    starting <- starting[-length(starting)] - 1L
    return(mapply("+", starting, lapply(jobsize, seq_len), SIMPLIFY=FALSE))
}

#################################################

.is_one_way <- function(design) 
# Checks if design matrix is a one-way layout.
{
    if (nrow(design) <= ncol(design)) {
        stop("design matrix has no residual degrees of freedom")
    }
    group <- designAsFactor(design)
    if (nlevels(group) == ncol(design)) {
        # Stripping out groups with only one level.
        groupings <- split(seq_len(nrow(design)), group)
        groupings[lengths(groupings)==1L] <- NULL
        return(groupings)
    } 
    return(NULL)
}

.ranksafe_qr <- function(design, tol=1e-7) 
# Rank-checking QR decomposition of a design matrix. Throws an
# error if the design matrix is not of full rank, which simplifies
# downstream processes as full rank can always be assumed.
{
    out <- qr(design, LAPACK=TRUE)
    d <- diag(out$qr)
    if (!all(abs(d) > tol)) { 
        stop("design matrix is not of full rank")
    }
    return(out)
}

#################################################

.calc_residuals_wt_zeroes <- function(x, design, QR, subset.row, lower.bound) 
# Computes residuals, but ensures that residuals for observations that 
# were below the lower bound are set to a constant value that is smaller 
# than all other residuals. Avoids spurious patterns when
# 
{
    if (!missing(design)) {
        QR <- .ranksafe_qr(design)
    }
    if (is.null(lower.bound)) { 
        stop("lower bound must be supplied or NA when computing residuals")
    }
    .Call(cxx_get_residuals, x, QR$qr, QR$qraux, subset.row - 1L, as.double(lower.bound))
}

.guess_lower_bound <- function(x, assay, lower.bound) 
# Getting the lower bound on the expression values for a given assay, if not supplied.
# We bump it up a little to make sure that expression values at the lower bound will 
# actually be detected as being "<= lower.bound".
{ 
    if (is.null(lower.bound)) { 
        if (assay=="exprs") {
            lower.bound <- log2(x@logExprsOffset) + 1e-8
        } else if (assay=="counts") {
            lower.bound <- 1e-8
        }
    }
    return(lower.bound)
}
