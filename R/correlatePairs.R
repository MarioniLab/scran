correlateNull <- function(ncells, iters=1e6, design=NULL, residuals=FALSE) 
# This builds a null distribution for the modified Spearman's rho.
#
# written by Aaron Lun
# created 10 February 2016
# last modified 27 May 2016
{
    if (!is.null(design)) { 
        if (!missing(ncells)) { 
            stop("cannot specify both 'ncells' and 'design'")
        }

        groupings <- .isOneWay(design)
        if (is.null(groupings) || residuals) { 
            # Using residualsd residual effects if the design matrix is not a one-way layout (or if forced by residuals=TRUE).
            QR <- .checkDesign(design)
            out <- .Call(cxx_get_null_rho_design, QR$qr, QR$qraux, as.integer(iters))
            if (is.character(out)) { 
                stop(out)
            }
        } else {
            # Otherwise, estimating the correlation as a weighted mean of the correlations in each group.
            # This avoids the need for the normality assumption in the residual effect simulation.
            out <- 0
            for (gr in groupings) {
                out.g <- .Call(cxx_get_null_rho, length(gr), as.integer(iters))
                if (is.character(out.g)) { 
                    stop(out.g)
                }
                out <- out + out.g * length(gr)
            }
            out <- out/nrow(design)
        }
        attrib <- list(design=design, residuals=residuals)

    } else {
        out <- .Call(cxx_get_null_rho, as.integer(ncells), as.integer(iters))
        if (is.character(out)) { 
            stop(out)
        }
        attrib <- NULL
    }

    # Storing attributes, to make sure it matches up.
    out <- sort(out)
    attributes(out) <- attrib
    return(out)  
}

.isOneWay <- function(design) {
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

setGeneric("correlatePairs", function(x, ...) standardGeneric("correlatePairs"))

setMethod("correlatePairs", "matrix", function(x, null.dist=NULL, design=NULL, BPPARAM=SerialParam(), use.names=TRUE, tol=1e-8, residuals=FALSE, subset.row=NULL)
# This calculates a (modified) Spearman's rho for each pair of genes.
#
# written by Aaron Lun
# created 10 February 2016
# last modified 6 August 2016
{
    compute.residuals <- FALSE
    if (!is.null(design)) { 
        QR <- .checkDesign(design)
        groupings <- .isOneWay(design)
        if (is.null(groupings) || residuals) { 
            compute.residuals <- TRUE
            groupings <- list(seq_len(ncol(x)))
        } 
        if (is.null(null.dist)) { 
            null.dist <- correlateNull(design=design, residuals=residuals)
        }
    } else {
        groupings <- list(seq_len(ncol(x)))
        if (is.null(null.dist)) { 
            null.dist <- correlateNull(ncol(x))
        } 
    }

    # Checking that the null distribution is sensible.
    if (!identical(design, attr(null.dist, "design"))) { 
        warning("'design' is not the same as that used to generate 'null.dist'")
    }
    if (!is.null(design)) { 
        if (!identical(residuals, attr(null.dist, "residuals"))) {
            warning("'residuals' is not the same as that used to generate 'null.dist'")
        }
    }
    null.dist <- as.double(null.dist)
    if (is.unsorted(null.dist)) { 
        null.dist <- sort(null.dist)
    }

    # Checking the subsetting and tolerance
    subset.row <- .subset_to_index(subset.row, x, byrow=TRUE)
    tol <- as.double(tol)

    # Generating all pairs of genes
    ngenes <- length(subset.row)
    if (ngenes < 2L) { stop("need at least two genes to compute correlations") }
    all.pairs <- combn(ngenes, 2L)
    gene1 <- all.pairs[1,]
    gene2 <- all.pairs[2,]

    # Iterating through all subgroups (for one-way layouts; otherwise, this is a loop of length 1).
    all.rho <- 0L
    for (subset.col in groupings) { 

        if (!compute.residuals) {
            # Ranking genes, in an error-tolerant way. This avoids getting untied rankings for zeroes
            # (which should have the same value +/- precision, as the prior count scaling cancels out).
            ranked.exprs <- .Call(cxx_rank_subset, x, subset.row - 1L, subset.col - 1L, tol)
        } else {
            # If we're computing residuals, we intervene here and replace values with the residuals.
            # Also need to replace the subset vector, as it'll already be subsetted.
            rx <- .Call(cxx_get_residuals, x, QR$qr, QR$qraux, subset.row - 1L)
            if (is.character(rx)) { stop(rx) }
            ranked.exprs <- .Call(cxx_rank_subset, rx, seq_len(nrow(rx)) - 1L, subset.col - 1L, tol)
        }
        if (is.character(ranked.exprs)) {
            stop(ranked.exprs)
        }

        # Running through each set of jobs 
        workass <- .workerAssign(length(gene1), BPPARAM)
        out <- bpmapply(FUN=.get_correlation, wstart=workass$start, wend=workass$end, BPPARAM=BPPARAM,
                        MoreArgs=list(gene1=gene1, gene2=gene2, ranked.exprs=ranked.exprs), SIMPLIFY=FALSE)
        current.rho <- unlist(out)

        # Adding a weighted value to the final.
        all.rho <- all.rho + current.rho * (length(subset.col)/ncol(x))
    }

    # Estimating the p-values (need to shift values to break ties conservatively by increasing the p-value).
    left <- findInterval(all.rho + 1e-8, null.dist)
    right <- length(null.dist) - findInterval(all.rho - 1e-8, null.dist)
    all.pval <- (pmin(left, right)+1)*2/(length(null.dist)+1)
    all.pval <- pmin(all.pval, 1)

    # Returning some useful output
    gene1 <- subset.row[gene1]
    gene2 <- subset.row[gene2]
    newnames <- NULL
    if (is.logical(use.names)) {
        if (use.names) {
            newnames <- rownames(x)
        }
    } else if (is.character(use.names)) {
        if (length(use.names)!=nrow(x)) {
            stop("length of 'use.names' does not match 'x' nrow")
        }
        newnames <- use.names
    }
    if (!is.null(newnames)) {
        gene1 <- newnames[gene1]
        gene2 <- newnames[gene2]
    }

    out <- data.frame(gene1=gene1, gene2=gene2, rho=all.rho, p.value=all.pval, 
                      FDR=p.adjust(all.pval, method="BH"), stringsAsFactors=FALSE)
    out <- out[order(out$p.value, -abs(out$rho)),]
    rownames(out) <- NULL
    return(out)
})

.workerAssign <- function(njobs, BPPARAM) {
    ncores <- bpworkers(BPPARAM)
    starting <- as.integer(seq(from=1, to=njobs+1, length.out=ncores+1))
    starting <- unique(starting[seq_len(ncores)])
    ending <- c((starting - 1L)[-1], njobs)
    return(list(start=starting, end=ending))
}

.get_correlation <- function(wstart, wend, gene1, gene2, ranked.exprs) {
    to.use <- wstart:wend
    out <- .Call(cxx_compute_rho, gene1[to.use], gene2[to.use], ranked.exprs)
    if (is.character(out)) { stop(out) }
    return(out)         
}

setMethod("correlatePairs", "SCESet", function(x, subset.row=NULL, ..., assay="exprs", get.spikes=FALSE) {
    if (is.null(subset.row)) {
        subset.row <- .spikeSubset(x, get.spikes)
    }
    correlatePairs(assayDataElement(x, assay), subset.row=subset.row, ...)             
})

