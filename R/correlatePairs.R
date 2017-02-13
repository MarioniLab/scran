setGeneric("correlatePairs", function(x, ...) standardGeneric("correlatePairs"))

.correlate_pairs <- function(x, null.dist=NULL, design=NULL, BPPARAM=SerialParam(), use.names=TRUE, tol=1e-8, 
                             iters=1e6, residuals=FALSE, subset.row=NULL, per.gene=FALSE)
# This calculates a (modified) Spearman's rho for each pair of genes.
#
# written by Aaron Lun
# created 10 February 2016
# last modified 2 February 2017
{
    compute.residuals <- FALSE
    if (!is.null(design)) { 
        QR <- qr(design, LAPACK=TRUE)
        groupings <- .isOneWay(design)
        if (is.null(groupings) || residuals) { 
            compute.residuals <- TRUE
            groupings <- list(seq_len(ncol(x)))
        } 
        if (is.null(null.dist)) { 
            null.dist <- correlateNull(design=design, residuals=residuals, iters=iters)
        }
    } else {
        groupings <- list(seq_len(ncol(x)))
        if (is.null(null.dist)) { 
            null.dist <- correlateNull(ncol(x), iters=iters)
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
    tol <- as.double(tol)
    pairings <- .construct_pair_indices(subset.row=subset.row, x=x)
    subset.row <- pairings$subset.row
    gene1 <- pairings$gene1
    gene2 <- pairings$gene2
    reorder <- pairings$reorder
    final.names <- .choose_gene_names(subset.row=subset.row, x=x, use.names=use.names)

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
                        MoreArgs=list(gene1=gene1 - 1L, gene2=gene2 - 1L, ranked.exprs=ranked.exprs), SIMPLIFY=FALSE)
        current.rho <- unlist(out)

        # Adding a weighted value to the final.
        all.rho <- all.rho + current.rho * (length(subset.col)/ncol(x))
    }

    # Estimating the p-values (need to shift values to break ties conservatively by increasing the p-value).
    left <- findInterval(all.rho + 1e-8, null.dist)
    right <- length(null.dist) - findInterval(all.rho - 1e-8, null.dist)
    limited <- left==0L | right==0L
    all.pval <- (pmin(left, right)+1)*2/(length(null.dist)+1)
    all.pval <- pmin(all.pval, 1)

    # Returning output on a per-gene basis, testing if each gene is correlated to any other gene.
    if (per.gene) {
        by.gene <- .Call(cxx_combine_corP, length(subset.row), gene1 - 1L, gene2 - 1L, 
                         all.rho, all.pval, limited, order(all.pval) - 1L) 
        if (is.character(by.gene)) stop(by.gene)
        out <- data.frame(gene=final.names, rho=by.gene[[2]], p.value=by.gene[[1]],
                          FDR=p.adjust(by.gene[[1]], method="BH"), 
                          limited=by.gene[[3]], stringsAsFactors=FALSE)
        rownames(out) <- NULL
        .is_sig_limited(out)
        return(out)
    }

    # Otherwise, returning the pairs themselves.
    gene1 <- final.names[gene1]
    gene2 <- final.names[gene2]
    out <- data.frame(gene1=gene1, gene2=gene2, rho=all.rho, p.value=all.pval, 
                      FDR=p.adjust(all.pval, method="BH"), 
                      limited=limited, stringsAsFactors=FALSE)
    if (reorder) {
        out <- out[order(out$p.value, -abs(out$rho)),]
        rownames(out) <- NULL
    }
    .is_sig_limited(out)
    return(out)
}

.construct_pair_indices <- function(subset.row, x) {
    is.ordered <- FALSE 
    if (is.matrix(subset.row)) {
        # If matrix, we're using pre-specified pairs.
        if (!is.numeric(subset.row) || ncol(subset.row)!=2L) { stop("'subset.row' should be a numeric matrix with 2 columns") }
        s1 <- .subset_to_index(subset.row[,1], x, byrow=TRUE)
        s2 <- .subset_to_index(subset.row[,2], x, byrow=TRUE)

        subset.row <- sort(unique(c(s1, s2)))
        m1 <- match(s1, subset.row)
        m2 <- match(s2, subset.row)
        gene1 <- pmin(m1, m2)
        gene2 <- pmax(m1, m2)
        is.ordered <- TRUE

    } else if (is.list(subset.row)) {
        # If list, we're correlating between one gene selected from each of two pools.
        if (length(subset.row)!=2L) { stop("'subset.row' as a list should have length 2") }
        converted <- lapply(subset.row, FUN=.subset_to_index, x=x, byrow=TRUE)
        converted <- lapply(converted, FUN=unique)
        if (any(lengths(converted)==0L)) { stop("need at least one gene in each set to compute correlations") }

        subset.row <- sort(unique(unlist(converted)))
        m1 <- match(converted[[1]], subset.row)
        m2 <- match(converted[[2]], subset.row)
        all.pairs <- expand.grid(m1, m2)
        gene1 <- pmin(all.pairs[,1], all.pairs[,2])
        gene2 <- pmax(all.pairs[,1], all.pairs[,2])

        # Pruning out redundant pairs. 
        if (length(intersect(m1, m2))!=0L) {
            o <- order(gene1, gene2)
            gene1 <- gene1[o]
            gene2 <- gene2[o]
            keep <- c(TRUE, diff(gene1)!=0L | diff(gene2)!=0L) & gene1!=gene2
            gene1 <- gene1[keep]
            gene2 <- gene2[keep]
        }

    } else {
        # Otherwise, it's assumed to be a single pool, and we're just correlating between pairs within it.
        subset.row <- .subset_to_index(subset.row, x, byrow=TRUE)
        ngenes <- length(subset.row)
        if (ngenes < 2L) { stop("need at least two genes to compute correlations") }
       
        # Generating all pairs of genes within the subset.
        all.pairs <- combn(ngenes, 2L)
        gene1 <- all.pairs[1,]
        gene2 <- all.pairs[2,]
    }

    return(list(subset.row=subset.row, gene1=gene1, gene2=gene2, reorder=!is.ordered))
}

.get_correlation <- function(wstart, wend, gene1, gene2, ranked.exprs) {
    to.use <- wstart:wend
    out <- .Call(cxx_compute_rho, gene1[to.use], gene2[to.use], ranked.exprs)
    if (is.character(out)) { stop(out) }
    return(out)         
}

.choose_gene_names <- function(subset.row, x, use.names) {
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
        subset.row <- newnames[subset.row]
    }
    return(subset.row)
}

.is_sig_limited <- function(results, threshold=0.05) {
    if (any(results$FDR > threshold & results$limited)) { 
        warning(sprintf("lower bound on p-values at a FDR of %s, increase 'iter'", as.character(threshold)))
    }
    invisible(NULL)
}

setMethod("correlatePairs", "matrix", .correlate_pairs)

setMethod("correlatePairs", "SCESet", function(x, subset.row=NULL, use.names=TRUE, per.gene=FALSE, ..., assay="exprs", get.spikes=FALSE) {
    by.spikes <- FALSE
    if (is.null(subset.row)) {
        subset.row <- .spikeSubset(x, get.spikes)
        by.spikes <- TRUE
    }
    out <- .correlate_pairs(assayDataElement(x, assay), subset.row=subset.row, per.gene=per.gene, use.names=use.names, ...)

    # Returning a row for all elements, even if it is NA.
    if (per.gene && by.spikes) {
        expanded <- rep(NA_integer_, nrow(x))
        expanded[subset.row] <- seq_len(nrow(out))
        out <- out[expanded,]
        out$gene <- .choose_gene_names(seq_len(nrow(x)), x, use.names)
        rownames(out) <- NULL
    }

    return(out)
})

