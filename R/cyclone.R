setGeneric("cyclone", function(x, ...) standardGeneric("cyclone"))

setMethod("cyclone", "ANY", function(x, pairs, gene.names=rownames(x), iter=1000, min.iter=100, min.pairs=50, BPPARAM=bpparam(), verbose=FALSE)
# Takes trained pairs and test data, and predicts the cell cycle phase from that. 
#
# written by Antonio Scialdone
# with modifications by Aaron Lun
# created 22 January 2016    
# last modified 17 February 2016
{ 
    x <- as.matrix(x)
    storage.mode(x) <- "double"
    Ngenes <- nrow(x)
    if (length(gene.names)!=Ngenes) {
        stop("length of 'gene.names' must be equal to 'x' nrows")
    }
    iter <- as.integer(iter)
    min.iter <- as.integer(min.iter)
    min.pairs <- as.integer(min.pairs)

    # Only keeping training pairs where both genes are in the test data.
    for (p in names(pairs)) {
        curp <- pairs[[p]]
        m1 <- match(curp$first, gene.names)
        m2 <- match(curp$second, gene.names)
        keep <- !is.na(m1) & !is.na(m2)
        m1 <- m1[keep]
        m2 <- m2[keep]
        pairs[[p]] <- list(first=m1, second=m2)
    }

    if (verbose) { 
        cat(sprintf("Number of G1 pairs: %d\n", nrow(pairs$G1)))
        cat(sprintf("Number of S pairs: %d\n", nrow(pairs$S)))
        cat(sprintf("Number of G2M pairs: %d\n", nrow(pairs$G2)))
    }
  
    # Run the allocation algorithm
    ncells <- ncol(x)
    workass <- .workerAssign(ncells, BPPARAM)
    all.cores <- seq_along(workass$start) 
    common.args <- list(X=all.cores, FUN=.get_phase_score, BPPARAM=BPPARAM,
                        work.start=workass$start, work.end=workass$end, 
                        Ngenes=Ngenes, exprs=x, iter=iter, min.iter=min.iter, min.pairs=min.pairs)
    out.G1 <- do.call(bplapply, c(list(pairings=pairs$G1), common.args))
    out.S <- do.call(bplapply, c(list(pairings=pairs$S), common.args))
    out.G2M <- do.call(bplapply, c(list(pairings=pairs$G2M), common.args))

    # Assembling the output.
    lapply(out.G1, FUN=function(y) { if (is.character(y)) { stop(y) } })
    lapply(out.S, FUN=function(y) { if (is.character(y)) { stop(y) } })
    lapply(out.G2M, FUN=function(y) { if (is.character(y)) { stop(y) } })
    score.G1 <- unlist(out.G1)
    score.S <- unlist(out.S)
    score.G2M <- unlist(out.G2M)
    
    scores <- data.frame(G1=score.G1, S=score.S, G2M=score.G2M)
    scores.normalised <- scores/rowSums(scores)
    return(list(scores=scores, normalized.scores=scores.normalised))  
})

.get_phase_score <- function(core, work.start, work.end, Ngenes, exprs, pairings, iter, min.iter, min.pairs) {
    to.use <- c(work.start[core], work.end[core])
    .Call(cxx_shuffle_scores, to.use, Ngenes, exprs, pairings$first, pairings$second, iter, min.iter, min.pairs) 
}

setMethod("cyclone", "SCESet", function(x, ..., assay="counts", get.spikes=FALSE) {
    cyclone(.getUsedMatrix(x, assay, get.spikes), ...)          
})

