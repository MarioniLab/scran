.overlapExprs <- function(x, groups, design=NULL, residuals=FALSE, tol=1e-8, subset.row=NULL)
# Computes the gene-specific overlap in expression profiles between two groups of cells.
# This aims to determine whether two distributions of expression values are well-separated.    
# 
# written by Aaron Lun
# created 17 April 2017
# last modified 18 April 2017
{
    compute.residuals <- FALSE
    if (!is.null(design)) { 
        QR <- qr(design, LAPACK=TRUE)
        groupings <- .isOneWay(design)
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
        use.x <- .Call(cxx_get_residuals, x, QR$qr, QR$qraux, subset.row - 1L)
        if (is.character(use.x)) { stop(use.x) }
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

    # Running through each blocking level and computing the proportions.
    for (subset.col in groupings) { 
        cur.groups <- groups[subset.col]
        by.group <- split(subset.col, cur.groups, drop=FALSE)
        if (length(by.group)==1L) { next }

        out <- .Call(cxx_overlap_exprs, use.x, use.subset.row, by.group, as.double(tol))
        if (is.character(out)) { stop(out) }
        output <- mapply("+", output, out[[1]], SIMPLIFY=FALSE)
        used.cells <- mapply("+", used.cells, out[[2]], SIMPLIFY=FALSE)
    }

    # Normalizing the output matrices.
    for (g in names(output)) { 
        output[[g]] <- t(t(output[[g]])/used.cells[[g]])
    }
    return(output)
}

setGeneric("overlapExprs", function(x, ...) standardGeneric("overlapExprs"))

setMethod("overlapExprs", "matrix", .overlapExprs)

setMethod("overlapExprs", "SCESet", function(x, ..., subset.row=NULL, assay="exprs", get.spikes=FALSE) {
    if (is.null(subset.row)) { subset.row <- .spikeSubset(x, get.spikes) }
    .overlapExprs(assayDataElement(x, assay), ..., subset.row=subset.row)
})                                 

