#' @importMethodsFrom BiocGenerics nrow ncol
#' @importFrom BiocGenerics colnames rownames
.check_batch_consistency <- function(batches, byrow=TRUE, ignore.null=FALSE) 
# Checking for identical number of rows (and rownames).
# It can also do the same for columns, if we're dealing with PC results.
{
    if (length(batches) < 2L) {
        return(NULL)
    }

    if (byrow) {
        DIMFUN <- nrow
        DIMNAMEFUN <- rownames
        DIM <- "row"
    } else {
        DIMFUN <- ncol
        DIMNAMEFUN <- colnames
        DIM <- "column"
    }

    first <- batches[[1]]
    ref.n <- DIMFUN(first)
    ref.names <- DIMNAMEFUN(first)

    for (b in 2:length(batches)) { 
        current <- batches[[b]]
        if (!identical(DIMFUN(current), ref.n)) {
            stop(sprintf("number of %ss is not the same across batches", DIM))
        }

        cur.names <- DIMNAMEFUN(current)
        if (ignore.null) { 
            if (is.null(cur.names)) { 
                cur.names <- ref.names
            } else if (is.null(ref.names)) {
                ref.names <- cur.names
            }
        }
        if (!identical(cur.names, ref.names)) {
            stop(sprintf("%s names are not the same across batches", DIM))
        }
    }

    return(NULL)
}

#' @importFrom SingleCellExperiment isSpike spikeNames
.check_spike_consistency <- function(batches) 
# Checking for identical spike-in sets and (overall) identities.
{
    if (length(batches) < 2L) {
        return(NULL)
    }

    ref.spike.names <- spikeNames(batches[[1]])
    ref.spike <- isSpike(batches[[1]])
    for (b in seq_along(batches)) {
        if (!identical(ref.spike.names, spikeNames(batches[[b]]))) {
            stop("spike-in sets differ across batches")
        }
        if (!identical(ref.spike, isSpike(batches[[b]]))) {
            stop("spike-in identities differ across batches")
        }
    }
    return(NULL)
}

#' @importFrom methods is
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
.SCEs_to_matrices <- function(batches, assay.type, subset.row, get.spikes)
# Convenience function to convert multiple SCEs to a list of expression matrices.
# Also redefines subset.row based on intersection with spike-ins, if necessary.    
{
    if (length(batches)<1L) {
        stop("at least one batch must be supplied")
    }
    .check_batch_consistency(batches, byrow=TRUE)

    all.sce <- vapply(batches, is, class2="SingleCellExperiment", FUN.VALUE=TRUE)
    if (length(unique(all.sce))!=1L) {
        stop("cannot mix SingleCellExperiments and other objects")
    }

    if (all(all.sce)) { 
        .check_spike_consistency(batches)
        subset.row <- .SCE_subset_genes(subset.row, batches[[1]], get.spikes)
        batches <- lapply(batches, assay, i=assay.type, withDimnames=FALSE)
    } else if (!is.null(subset.row)) {
        subset.row <- .subset_to_index(subset.row, batches[[1]], byrow=TRUE)
    }
	return(list(batches=batches, subset.row=subset.row))
}
