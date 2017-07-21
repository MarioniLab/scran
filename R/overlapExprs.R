.overlapExprs <- function(x, groups, design=NULL, residuals=FALSE, tol=1e-8, 
                          BPPARAM=SerialParam(), subset.row=NULL, lower.bound=NULL)
# Computes the gene-specific overlap in expression profiles between two groups of cells.
# This aims to determine whether two distributions of expression values are well-separated.    
# 
# written by Aaron Lun
# created 17 April 2017
# last modified 27 April 2017
{
    compute.residuals <- FALSE
    if (!is.null(design)) { 
        groupings <- .is_one_way(design)
        if (is.null(groupings) || residuals) { 
            compute.residuals <- TRUE
            groupings <- list(seq_len(ncol(x)))
        } 
    } else {
        groupings <- list(seq_len(ncol(x)))
    }

    # Checking dimensions.
    subset.row <- .subset_to_index(subset.row, x, byrow=TRUE)
    ncells <- ncol(x)
    if (length(groups)!=ncells) { 
        stop("length of 'groups' not equal to number of cells")
    }
    if (!is.null(design) && nrow(design)!=ncells) { 
        stop("'nrow(design)' not equal to number of cells")
    }

    # Computing residuals; also replacing the subset vector, as it'll already be subsetted.
    if (compute.residuals) { 
        use.x <- .calc_residuals_wt_zeroes(x, design, subset.row=subset.row, lower.bound=lower.bound)
        use.subset.row <- seq_len(nrow(use.x)) - 1L
    } else {
        use.x <- x
        use.subset.row <- subset.row - 1L
    }

    # Setting up the output matrices.
    groups <- as.factor(groups)
    unique.groups <- levels(groups)
    ngroups <- nlevels(groups)
    output <- used.cells <- vector("list", ngroups)
    ngenes <- length(subset.row)
    for (g in seq_len(ngroups)) { 
        temp <- matrix(0, ngenes, ngroups-1L) 
        colnames(temp) <- unique.groups[-g]
        rownames(temp) <- rownames(x)[subset.row]
        output[[g]] <- temp
        temp.n <- integer(ngroups-1L)
        names(temp.n) <- colnames(temp)
        used.cells[[g]] <- temp.n
    }
    names(output) <- names(used.cells) <- unique.groups

    # Split up the jobs by genes for multicore execution
    wout <- .worker_assign(length(use.subset.row), BPPARAM)
    for (i in seq_along(wout)) {
        wout[[i]] <- use.subset.row[wout[[i]]]
    }
    tol <- as.double(tol)

    # Running through each blocking level, splitting cells into groups and computing the proportions.
    for (subset.col in groupings) { 
        cur.groups <- groups[subset.col]
        by.group <- split(subset.col, cur.groups, drop=FALSE)
        if (length(by.group)==1L) { next }

        bout <- bplapply(wout, FUN=.find_overlap_exprs, x=use.x, by.group=by.group, tol=tol, BPPARAM=BPPARAM)

        # Adding the proportions to what we already have. 
        props <- do.call(rbind, lapply(bout, "[[", i=1)) 
        numsc <- bout[[1]][[2]]
        output <- mapply("+", output, props, SIMPLIFY=FALSE)
        used.cells <- mapply("+", used.cells, numsc, SIMPLIFY=FALSE)
    }

    # Normalizing the output matrices.
    for (g in names(output)) { 
        output[[g]] <- t(t(output[[g]])/used.cells[[g]])
    }
    return(output)
}

.find_overlap_exprs <- function(x, subset.row, by.group, tol) 
# Pass all arguments explicitly rather than via function environment
# (avoid duplication of memory in bplapply).
{
    .Call(cxx_overlap_exprs, x, subset.row, by.group, tol)
}

setGeneric("overlapExprs", function(x, ...) standardGeneric("overlapExprs"))

setMethod("overlapExprs", "matrix", .overlapExprs)

setMethod("overlapExprs", "SCESet", function(x, ..., subset.row=NULL, lower.bound=NULL, assay="exprs", get.spikes=FALSE) {
    if (is.null(subset.row)) { subset.row <- .spike_subset(x, get.spikes) }
    lower.bound <- .guess_lower_bound(x, assay, lower.bound)
    .overlapExprs(assayDataElement(x, assay), ..., lower.bound=lower.bound, subset.row=subset.row)
})                                 

