.overlapExprs <- function(x, groups, 
                          block=NULL, 
                          design=NULL, 
                          rank.type=c("any", "all"),
                          direction=c("any", "up", "down"),
                          tol=1e-8, 
                          BPPARAM=SerialParam(), 
                          subset.row=NULL, 
                          residuals=FALSE, 
                          lower.bound=NULL)
# Computes the gene-specific overlap in expression profiles between two groups of cells.
# This aims to determine whether two distributions of expression values are well-separated.    
# 
# written by Aaron Lun
# created 17 April 2017
{
    subset.row <- .subset_to_index(subset.row, x, byrow=TRUE)
    tol <- as.double(tol)
    ncells <- ncol(x)
    if (length(groups)!=ncells) { 
        stop("length of 'groups' not equal to number of cells")
    }
    if (!is.null(design) && nrow(design)!=ncells) { 
        stop("'nrow(design)' not equal to number of cells")
    }

    # Computing residuals; also replacing the subset vector, as it'll already be subsetted.
    by.block <- .check_blocking_inputs(x, block=block, design=design, residuals=residuals)
    if (!is.null(design) && is.null(block)) { 
        use.x <- .calc_residuals_wt_zeroes(x, design, subset.row=subset.row, lower.bound=lower.bound)
        use.subset.row <- seq_len(nrow(use.x)) - 1L
    } else {
        use.x <- x
        use.subset.row <- subset.row - 1L
    }

    # Setting up the output matrices.
    groups <- as.factor(groups)
    overlap <- used.cells <- rep(list(0), nlevels(groups))
    names(overlap) <- names(used.cells) <- levels(groups) 

    # Split up the jobs by genes for multicore execution
    wout <- .worker_assign(length(use.subset.row), BPPARAM)
    for (i in seq_along(wout)) {
        wout[[i]] <- use.subset.row[wout[[i]]]
    }

    # Running through each blocking level, splitting cells into groups.
    for (subset.col in by.block) { 
        cur.groups <- groups[subset.col]
        by.group <- split(subset.col, cur.groups, drop=FALSE)
        stopifnot(identical(names(by.group), levels(groups))) # Sanity check.

        # Computing the proportions within this block.
        bout <- bplapply(wout, FUN=.find_overlap_exprs, x=use.x, by.group=by.group, tol=tol, BPPARAM=BPPARAM)

        # Adding the proportions to what we already have. 
        props <- do.call(rbind, lapply(bout, "[[", i=1)) 
        numsc <- bout[[1]][[2]]
        overlap <- mapply("+", overlap, props, SIMPLIFY=FALSE)
        used.cells <- mapply("+", used.cells, numsc, SIMPLIFY=FALSE)
    }

    # Normalizing the overlap matrices.
    for (g in names(overlap)) { 
        overlap[[g]] <- t(t(overlap[[g]])/used.cells[[g]]) 
    }

    rank.type <- match.arg(rank.type)
    direction <- match.arg(direction)
    gene.names <- rownames(x)[subset.row]
    if (is.null(gene.names)) { 
        gene.names <- subset.row # using indices instead.
    }
    .prepare_overlap_output(overlap, rank.type=rank.type, direction=direction, gene.names=gene.names)
}

###########################################################
# Internal functions.
###########################################################

.check_blocking_inputs <- function(x, block, design, residuals) 
# This makes sure that the blocking inputs are sensible.
{
    if (!is.null(block)) { 
        blocks <- split(seq_len(ncol(x)), block)
    } else { 
        blocks <- list(seq_len(ncol(x)))

        # And a series of warnings...
        if (!is.null(design)) { 
            if (is.null(.is_one_way(design))) { 
                if (residuals) {
                    .Deprecated(msg="'residuals=TRUE' is deprecated, choose between 'design' and 'block'")
                }
            } else if (!residuals) {
                .Deprecated(msg="'residuals=FALSE' is deprecated, use 'block' instead")
            }
        }
    }
    return(blocks)
}

.find_overlap_exprs <- function(x, subset.row, by.group, tol) 
# Pass all arguments explicitly rather than via function environment
# (avoid duplication of memory in bplapply).
{
    .Call(cxx_overlap_exprs, x, subset.row, by.group, tol)
}

.prepare_overlap_output <- function(overlaps, direction, rank.type, gene.names) 
# This handles the final part of the function; ordering and 
# preparing the output values.
{
    output <- vector("list", length(overlaps))
    names(output) <- names(overlaps)
    for (i in seq_along(overlaps)) { 
        current <- overlaps[[i]]

        # Figuring out the direction to use for ranking.
        if (direction=="any") {
            metric <- 0.5-abs(current-0.5)
        } else if (direction=="up") {
            metric <- 1-current
        } else {
            metric <- current
        }
        if (!ncol(metric)) { 
            metric <- matrix(NA_real_, nrow(metric), 1L)
        }
       
        # Figuring out how exactly to rank them. 
        if (rank.type=="any") {
            rank.out <- .rank_top_genes(metric)
            to.add <- rank.out$rank
            field.name <- "Top"
            o <- order(to.add, rank.out$value)
        } else {
            # Getting the smallest overlap and ranking by that.
            to.add <- current[.find_largest_index(metric)]
            field.name <- "Worst"
            o <- order(to.add)
        }

        # Saving everything as a data frame.
        odata <- DataFrame(to.add, current, row.names=gene.names)
        colnames(odata) <- c(field.name, sprintf("overlap.%s", names(overlaps)[-i]))
        output[[i]] <- odata[o,,drop=FALSE]
    }
    return(output)
}

###########################################################
# S4 method definitions
###########################################################

setGeneric("overlapExprs", function(x, ...) standardGeneric("overlapExprs"))

setMethod("overlapExprs", "ANY", .overlapExprs)

setMethod("overlapExprs", "SingleCellExperiment", 
          function(x, ..., subset.row=NULL, lower.bound=NULL, assay.type="logcounts", get.spikes=FALSE) {

    subset.row <- .SCE_subset_genes(subset.row=subset.row, x=x, get.spikes=get.spikes)
    lower.bound <- .guess_lower_bound(x, assay.type, lower.bound)
    .overlapExprs(assay(x, i=assay.type), ..., lower.bound=lower.bound, subset.row=subset.row)
})                                 

