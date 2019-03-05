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

    return(list(is.spike=is.spike, sf.cell=sf.cell, sf.spike=sf.spike))
} 
