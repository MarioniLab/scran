#' @importFrom BiocParallel bplapply SerialParam
.correlate_pairs <- function(x, null.dist=NULL, tol=1e-8, iters=1e6, 
                             block=NULL, design=NULL, residuals=FALSE, lower.bound=NULL, 
                             use.names=TRUE, subset.row=NULL, pairings=NULL, per.gene=FALSE, 
                             cache.size=100L, BPPARAM=SerialParam())
# This calculates a (modified) Spearman's rho for each pair of genes.
#
# written by Aaron Lun
# created 10 February 2016
{
    null.out <- .check_null_dist(x, block=block, design=design, iters=iters, null.dist=null.dist, residuals=residuals)
    null.dist <- null.out$null
    by.block <- null.out$blocks

    # Checking which pairwise correlations should be computed.
    cache.size <- as.integer(cache.size)
    pair.out <- .construct_pair_indices(subset.row=subset.row, x=x, pairings=pairings)
    subset.row <- pair.out$subset.row
    gene1 <- pair.out$gene1
    gene2 <- pair.out$gene2
    reorder <- pair.out$reorder

    # Computing residuals (setting values that were originally zero to a lower bound).
    # Also replacing the subset vector, as it'll already be subsetted.
    if (!is.null(design) && is.null(block)) {
        use.x <- .calc_residuals_wt_zeroes(x, design, subset.row=subset.row, lower.bound=lower.bound) 
        use.subset.row <- seq_len(nrow(use.x)) - 1L
    } else {
        use.x <- x
        use.subset.row <- subset.row - 1L
    }

    # Splitting up gene pairs into jobs for multicore execution, converting to 0-based indices.
    wout <- .worker_assign(length(gene1), BPPARAM)
    sgene1 <- sgene2 <- vector("list", length(wout))
    for (i in seq_along(wout)) {
        sgene1[[i]] <- gene1[wout[[i]]] - 1L 
        sgene2[[i]] <- gene2[wout[[i]]] - 1L
    }

    # Iterating through all blocking levels (for one-way layouts; otherwise, this is a loop of length 1).
    # Computing correlations between gene pairs, and adding a weighted value to the final average.
    all.rho <- numeric(length(gene1))
    for (subset.col in by.block) { 
        ranked.exprs <- .Call(cxx_rank_subset, use.x, use.subset.row, subset.col - 1L, as.double(tol))
        out <- bpmapply(FUN=.get_correlation, gene1=sgene1, gene2=sgene2, 
                        MoreArgs=list(ranked.exprs=ranked.exprs, cache.size=cache.size), 
                        BPPARAM=BPPARAM, SIMPLIFY=FALSE)
        current.rho <- unlist(out)
        all.rho <- all.rho + current.rho * length(subset.col)/ncol(x)
    }

    # Estimating the p-values (need to shift values to break ties conservatively by increasing the p-value).
    left <- findInterval(all.rho + 1e-8, null.dist)
    right <- length(null.dist) - findInterval(all.rho - 1e-8, null.dist)
    limited <- left==0L | right==0L
    all.pval <- (pmin(left, right)+1)*2/(length(null.dist)+1)
    all.pval <- pmin(all.pval, 1)

    # Returning output on a per-gene basis, testing if each gene is correlated to any other gene.
    final.names <- .choose_gene_names(subset.row=subset.row, x=x, use.names=use.names)
    if (per.gene) {
        by.gene <- .Call(cxx_combine_corP, length(subset.row), gene1 - 1L, gene2 - 1L, 
                         all.rho, all.pval, limited, order(all.pval) - 1L) 

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
    out <- DataFrame(gene1=gene1, gene2=gene2, rho=all.rho, p.value=all.pval, 
                     FDR=p.adjust(all.pval, method="BH"), 
                     limited=limited, stringsAsFactors=FALSE)
    if (reorder) {
        out <- out[order(out$p.value, -abs(out$rho)),]
        rownames(out) <- NULL
    }
    .is_sig_limited(out)
    return(out)
}

.check_null_dist <- function(x, block, design, iters, null.dist, residuals) 
# This makes sure that the null distribution is in order.
{
    if (!is.null(block)) { 
        blocks <- split(seq_len(ncol(x)), block)
        if (is.null(null.dist)) { 
            null.dist <- correlateNull(block=block, iters=iters)
        }

    } else if (!is.null(design)) { 
        if (is.null(.is_one_way(design))) { 
            if (residuals) {
                .Deprecated(msg="'residuals=TRUE' is deprecated, choose between 'design' and 'block'")
            }
        } else if (!residuals) {
            .Deprecated(msg="'residuals=FALSE' is deprecated, use 'block' instead")
        }
            
        blocks <- list(seq_len(ncol(x)))
        if (is.null(null.dist)) { 
            null.dist <- correlateNull(design=design, iters=iters)
        }

    } else {
        blocks <- list(seq_len(ncol(x)))
        if (is.null(null.dist)) { 
            null.dist <- correlateNull(ncol(x), iters=iters)
        } 
    }

    # Checking that the null distribution is sensible.
    if (!identical(block, attr(null.dist, "block"))) { 
        warning("'block' is not the same as that used to generate 'null.dist'")
    }
    if (!identical(design, attr(null.dist, "design"))) { 
        warning("'design' is not the same as that used to generate 'null.dist'")
    }
    null.dist <- as.double(null.dist)
    if (is.unsorted(null.dist)) { 
        null.dist <- sort(null.dist)
    }
    
    return(list(null=null.dist, blocks=blocks))
}

#' @importFrom utils combn
.construct_pair_indices <- function(subset.row, x, pairings) 
# This returns a new subset-by-row vector, along with the pairs of elements
# indexed along that vector (i.e., "1" refers to the first element of subset.row,
# rather than the first element of "x").
{
    subset.row <- .subset_to_index(subset.row, x, byrow=TRUE)
    reorder <- TRUE

    if (is.matrix(pairings)) {
        # If matrix, we're using pre-specified pairs.
        if ((!is.numeric(pairings) && !is.character(pairings)) || ncol(pairings)!=2L) { 
            stop("'pairings' should be a numeric/character matrix with 2 columns") 
        }
        s1 <- .subset_to_index(pairings[,1], x, byrow=TRUE)
        s2 <- .subset_to_index(pairings[,2], x, byrow=TRUE)

        # Discarding elements not in subset.row.
        keep <- s1 %in% subset.row & s2 %in% subset.row
        s1 <- s1[keep]
        s2 <- s2[keep]

        subset.row <- sort(unique(c(s1, s2)))
        m1 <- match(s1, subset.row)
        m2 <- match(s2, subset.row)
        gene1 <- pmin(m1, m2)
        gene2 <- pmax(m1, m2)
        reorder <- FALSE

    } else if (is.list(pairings)) {
        # If list, we're correlating between one gene selected from each of two pools.
        if (length(pairings)!=2L) { 
            stop("'pairings' as a list should have length 2") 
        }
        converted <- lapply(pairings, FUN=function(gene.set) {
            gene.set <- .subset_to_index(gene.set, x=x, byrow=TRUE)
            intersect(gene.set, subset.row) # automatically gets rid of duplicates.
        })
        if (any(lengths(converted)==0L)) { 
            stop("need at least one gene in each set to compute correlations") 
        }

        subset.row <- sort(unique(unlist(converted)))
        m1 <- match(converted[[1]], subset.row)
        m2 <- match(converted[[2]], subset.row)
        all.pairs <- expand.grid(m1, m2)
        gene1 <- all.pairs[,1]
        gene2 <- all.pairs[,2]

    } else if (is.null(pairings)) {
        # Otherwise, it's assumed to be a single pool, and we're just correlating between pairs within it.
        ngenes <- length(subset.row)
        if (ngenes < 2L) { 
            stop("need at least two genes to compute correlations") 
        }
       
        # Generating all pairs of genes within the subset.
        all.pairs <- combn(ngenes, 2L)
        gene1 <- all.pairs[1,]
        gene2 <- all.pairs[2,]

    } else {
        stop("pairings should be a list, matrix or NULL")
    }

    return(list(subset.row=subset.row, gene1=gene1, gene2=gene2, reorder=reorder))
}

.get_correlation <- function(gene1, gene2, ranked.exprs, cache.size) 
# Pass all arguments explicitly rather than through the function environments
# (avoid duplicating memory in bplapply).
{
    .Call(cxx_compute_rho, gene1, gene2, ranked.exprs, cache.size)
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

#' @export
setGeneric("correlatePairs", function(x, ...) standardGeneric("correlatePairs"))

#' @export
setMethod("correlatePairs", "ANY", .correlate_pairs)

#' @importFrom SummarizedExperiment assay
#' @export
setMethod("correlatePairs", "SingleCellExperiment", 
          function(x, ..., use.names=TRUE, subset.row=NULL, per.gene=FALSE, 
                   lower.bound=NULL, assay.type="logcounts", get.spikes=FALSE) {

    subset.row <- .SCE_subset_genes(subset.row, x=x, get.spikes=get.spikes)              
    lower.bound <- .guess_lower_bound(x, assay.type, lower.bound)
    .correlate_pairs(assay(x, i=assay.type), subset.row=subset.row, per.gene=per.gene, 
                     use.names=use.names, lower.bound=lower.bound, ...)
})

