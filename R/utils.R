#' @importFrom SingleCellExperiment isSpike
.spike_subset <- function(x, get.spikes) {
    if (!get.spikes) {
        nokeep <- isSpike(x)
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

.SCE_subset_genes <- function(subset.row, x, get.spikes) {
    despiked <- .spike_subset(x, get.spikes)
    if (is.null(subset.row)) { 
        subset.row <- despiked
    } else {
        subset.row <- .subset_to_index(subset.row, x, byrow=TRUE)
        if (!is.null(despiked)) { 
            subset.row <- intersect(subset.row, which(despiked))
        }
    }
    return(subset.row)
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
    if (nrow(design)!=ncol(x)) {
        stop("number of rows in 'design' should be equal to 'ncol(x)'")
    }
    return(list(design=design))
}

#' @importFrom SummarizedExperiment assayNames
#' @importFrom scater librarySizeFactors
#' @importFrom SingleCellExperiment sizeFactorNames
#' @importFrom BiocGenerics sizeFactors
.check_centered_SF <- function(x, assay.type, block=NULL)
# Checks if 'logcounts' was requested, and if it could have been computed from counts.
# If so, then it checks whether the size factors are centered across gene sets.
{
    if (assay.type=="logcounts" && "counts" %in% assayNames(x)) {
        if (is.null(block)) {
            by.block <- list(seq_len(ncol(x)))
        } else {
            by.block <- split(seq_len(ncol(x)), block)
        }

        FUN <- function(sfs) {
            means <- numeric(length(by.block))
            for (idx in seq_along(by.block)) {
                means[idx] <- mean(sfs[by.block[[idx]]])
            }
            return(means)
        }
        
        ref <- sizeFactors(x)
        if (is.null(ref)) {
            ref <- librarySizeFactors(x)
        }
        ref.means <- FUN(ref)

        is.okay <- TRUE
        for (sf in sizeFactorNames(x)) {
            sf.means <- FUN(sizeFactors(x, sf))
            if (!isTRUE(all.equal(ref.means, sf.means))) {
                is.okay <- FALSE
                break
            }
        }

        if (!is.okay) {
            if (is.null(block)) { 
                warning("size factors not centred, run 'normalize()' first")
            } else {
                warning("size factors not centred, run 'multiBlockNorm()' first")
            }
        }
    }
    return(NULL)
}

#' @importFrom BiocGenerics sizeFactors
#' @importFrom SingleCellExperiment spikeNames isSpike
.prepare_cv2_data <- function(x, spike.type) 
# Prepares data for calculation of CV2.
# In particular, extracting spike-ins and their size factors.    
{
    sf.cell <- sizeFactors(x)
    if (is.null(spike.type) || !is.na(spike.type)) { 
        if (is.null(spike.type)) { 
            # Get all spikes.
            spike.type <- spikeNames(x)            
        } else {
            if (!all(spike.type %in% spikeNames(x))) { 
                stop(sprintf("spike-in set '%s' does not exist", spike.type[1]))
            }
        }
        if (!length(spike.type)) { 
            stop("no spike-in sets specified from 'x'")
        }

        # Collecting the size factors for the requested spike-in sets.
        # Check that all spike-in factors are either NULL or identical.
        collected <- NULL
        is.spike <- logical(nrow(x))
        for (st in seq_along(spike.type)) {
            cur.type <- spike.type[st]
            is.spike <- is.spike | isSpike(x, type=cur.type)
            cur.sf <- sizeFactors(x, type=cur.type)
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

#' @importFrom BiocParallel bpnworkers
.worker_assign <- function(njobs, BPPARAM)
# Assigns jobs to workers, where the each element of the output list is 
# the index vector of the jobs. These are guaranteed to be consecutive
# so the bplapply's output can just be directly combined.
{
    ncores <- bpnworkers(BPPARAM)
    starting <- as.integer(seq(from=1, to=njobs+1, length.out=ncores+1))
    jobsize <- diff(starting)
    starting <- starting[-length(starting)] - 1L
    return(mapply("+", starting, lapply(jobsize, seq_len), SIMPLIFY=FALSE))
}

.split_vector_by_workers <- function(vec, assignments) 
# Convenience function to split a vector by the assigned core,
# to get something that can be easily used in 'bplapply'.
{
    by.core <- vector("list", length(assignments))
    for (core in seq_along(assignments)) {
        by.core[[core]] <- vec[assignments[[core]]]
    }
    names(by.core) <- names(assignments)
    return(by.core)
}

.split_matrix_by_workers <- function(mat, assignments, byrow=TRUE) 
# Convenience function to split a mat by the assigned core,
# to get something that can be easily used in 'bplapply'.
{
    by.core <- vector("list", length(assignments))
    for (core in seq_along(assignments)) {
        chosen <- assignments[[core]]
        if (byrow) {
            current <- mat[chosen,,drop=FALSE]
        } else {
            current <- mat[,chosen,drop=FALSE]
        }
        by.core[[core]] <- current
    }
    names(by.core) <- names(assignments)
    return(by.core)
}

#################################################

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

.guess_lower_bound <- function(x, assay.type, lower.bound) 
# Getting the lower bound on the expression values for a given assay, if not supplied.
# We bump it up a little to make sure that expression values at the lower bound will 
# actually be detected as being "<= lower.bound".
{ 
    if (is.null(lower.bound)) { 
        if (assay.type=="logcounts") {
            lower.bound <- log2(.get_log_offset(x)) + 1e-8
        } else if (assay.type=="counts") {
            lower.bound <- 1e-8
        }
    }
    return(lower.bound)
}

#' @importFrom S4Vectors metadata
.get_log_offset <- function(x) 
# Helper function to get the log-offset value.
{
    metadata(x)$log.exprs.offset
}
