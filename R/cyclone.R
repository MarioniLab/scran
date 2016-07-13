setGeneric("cyclone", function(x, ...) standardGeneric("cyclone"))

setMethod("cyclone", "matrix", function(x, pairs, gene.names=rownames(x), iter=1000, min.iter=100, min.pairs=50, BPPARAM=bpparam(), verbose=FALSE, subset.row=NULL)
# Takes trained pairs and test data, and predicts the cell cycle phase from that. 
#
# written by Antonio Scialdone
# with modifications by Aaron Lun
# created 22 January 2016    
# last modified 7 June 2016
{ 
    if (length(gene.names)!=nrow(x)) {
        stop("length of 'gene.names' must be equal to 'x' nrows")
    }
    iter <- as.integer(iter)
    min.iter <- as.integer(min.iter)
    min.pairs <- as.integer(min.pairs)
   
    # Checking subset vector and blanking out the unused names.
    subset.row <- .subset_to_index(subset.row, x, byrow=TRUE)
    gene.names[-subset.row] <- NA
    
    # Only keeping training pairs where both genes are in the test data.
    for (p in names(pairs)) {
        curp <- pairs[[p]]
        m1 <- match(curp$first, gene.names)
        m2 <- match(curp$second, gene.names)
        keep <- !is.na(m1) & !is.na(m2)
        m1 <- m1[keep]
        m2 <- m2[keep]

        # Reformatting it to be a bit easier to access during permutations.
        retained <- logical(length(gene.names))
        retained[m1] <- TRUE
        retained[m2] <- TRUE
        new.indices <- cumsum(retained)
        pairs[[p]] <- list(first=new.indices[m1]-1L, second=new.indices[m2]-1L, 
                           index=which(retained)-1L) # For zero indexing.
    }

    if (verbose) { 
        cat(sprintf("Number of G1 pairs: %d\n", length(pairs$G1[[1]])))
        cat(sprintf("Number of S pairs: %d\n", length(pairs$S[[1]])))
        cat(sprintf("Number of G2M pairs: %d\n", length(pairs$G2[[1]])))
    }
  
    # Run the allocation algorithm
    ncells <- ncol(x)
    workass <- .workerAssign(ncells, BPPARAM)
    all.cores <- seq_along(workass$start) 
    common.args <- list(X=all.cores, FUN=.get_phase_score, BPPARAM=BPPARAM,
                        work.start=workass$start, work.end=workass$end, 
                        exprs=x, iter=iter, min.iter=min.iter, min.pairs=min.pairs)
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

    # Getting the phases.
    phases <- ifelse(score.G1 >= score.G2M, "G1", "G2M")
    phases[score.G1 < 0.5 & score.G2M < 0.5] <- "S"

    return(list(phases=phases, scores=scores, normalized.scores=scores.normalised))  
})

.get_phase_score <- function(core, work.start, work.end, exprs, pairings, iter, min.iter, min.pairs) {
    to.use <- c(work.start[core], work.end[core])
    .Call(cxx_shuffle_scores, to.use, exprs, pairings$first, pairings$second, pairings$index, iter, min.iter, min.pairs) 
}

setMethod("cyclone", "SCESet", function(x, pairs, subset.row=NULL, ..., assay="counts", get.spikes=FALSE) {
    if (is.null(subset.row)) {
        subset.row <- .spikeSubset(x, get.spikes)
    }
    cyclone(assayDataElement(x, assay), pairs=pairs, subset.row=subset.row, ...)          
})

