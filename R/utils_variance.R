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
    if (is.null(spike.type) || all(!is.na(spike.type))) { 
        if (is.null(spike.type)) { 
            # Get all spikes.
            spike.type <- spikeNames(x)            
        } else if (!all(spike.type %in% spikeNames(x))) { 
            stop(sprintf("spike-in set '%s' does not exist", spike.type[1]))
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

    list(is.spike=is.spike, sf.cell=sf.cell, sf.spike=sf.spike)
}

#' @importFrom BiocParallel bplapply
.lognormvar <- function(x, block, subset.row, design, BPPARAM, size.factors, pseudo.count) {
    subset.row <- .subset_to_index(subset.row, x, byrow=TRUE)
    wout <- .worker_assign(length(subset.row), BPPARAM)
    by.core <- .split_vector_by_workers(subset.row, wout)
    by.core <- .split_matrix_by_workers(x, by.core)

    if (!is.null(block)) { 
        if (ncol(x)!=length(block)) {
            stop("length of 'block' should be the same as 'ncol(x)'")
        }

        # Checking residual d.f.
        by.block <- split(seq_len(ncol(x))-1L, block, drop=TRUE)
	} else {
        by.block <- list(seq_len(ncol(x))-1L)
	}

    resid.df <- lengths(by.block) - 1L
	if (all(resid.df<=0L)){ 
		stop("no residual d.f. in any level of 'block' for variance estimation")
	}

    if (is.null(design)) {
	    raw.stats <- bplapply(by.core, FUN=compute_blocked_stats_lognorm, bygroup=by.block, 
            sf=size.factors, pseudo=pseudo.count, BPPARAM=BPPARAM)
        means <- do.call(rbind, lapply(raw.stats, FUN=function(x) t(x[[1]])))
        vars <- do.call(rbind, lapply(raw.stats, FUN=function(x) t(x[[2]])))
    } else {
        # Put linear modelling section here.
    }

	dimnames(means) <- dimnames(vars) <- list(rownames(x)[subset.row], names(by.block))

    list(means=means, vars=vars)
}

#' @importFrom BiocParallel SerialParam
#' @importFrom scater librarySizeFactors
.compute_var_stats_from_counts <- function(x, size.factors=NULL, subset.row=NULL, block=NULL, 
    fit.x=NULL, fit.size.factors=NULL, BPPARAM=SerialParam(), FUN, ...) 
{
    if (is.null(size.factors)) {
        size.factors <- librarySizeFactors(x, subset_row=subset.row)        
    }
    stats.out <- FUN(x, block=block, ..., subset.row=subset.row, BPPARAM=BPPARAM, size.factors=size.factors)

    if (is.null(fit.x)) {
        fit.stats <- stats.out
    } else {
        if (is.null(fit.size.factors)) {
            fit.size.factors <- librarySizeFactors(fit.x) # no subset_row here, as that only applies to 'x'.
        }

        # Rescaling so that the mean fit.size.factors is the same as each size.factors in each block.
        if (is.null(block)) {
            fit.size.factors <- fit.size.factors/mean(fit.size.factors) * mean(size.factors)
        } else {
            by.block <- split(seq_along(block), block)
            for (i in by.block) {
                current <- fit.size.factors[i]
                fit.size.factors[i] <- current / mean(current) * mean(size.factors[i])
            }
        }

        fit.stats <- FUN(fit.x, block=block, design=design, subset.row=subset.row, 
            BPPARAM=BPPARAM, size.factors=fit.size.factors)
    }

    list(x=stats.out, fit=fit.stats)
}

#' @importFrom stats density approx
.inverse_density_weights <- function(x, adjust=1) {
    out <- density(x, adjust=adjust, from=min(x), to=max(x))
    w <- 1/approx(out$x, out$y, xout=x)$y 
    w/mean(w)
}

#' @importFrom limma weighted.median
.correct_logged_expectation <- function(x, y, w, FUN) 
# Adjusting for any scale shift due to fitting to the log-values.
# The expectation of the log-values should be the log-expectation
# plus a factor that is dependent on the variance of the raw values
# divided by the squared mean, using a second-order Taylor approximation. 
# If we assume that the standard deviation of the variances is proportional
# to the mean variances with the same constant across all abundances,
# we should be able to correct the discrepancy with a single rescaling factor. 
{
    leftovers <- y/FUN(x)
    med <- weighted.median(leftovers, w, na.rm=TRUE)

    OUT <- function(x) { 
        output <- FUN(x) * med
        names(output) <- names(x)
        output
    }

    # We assume ratios are normally distributed around 1 with some standard deviation.
    std.dev <- unname(weighted.median(abs(leftovers/med - 1), w, na.rm=TRUE)) * 1.4826 
    list(trend=OUT, std.dev=std.dev)
}
