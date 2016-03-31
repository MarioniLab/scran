correlateNull <- function(ncells, iters=1e6, design=NULL) 
# This builds a null distribution for the modified Spearman's rho.
#
# written by Aaron Lun
# created 10 February 2016
# last modified 21 February 2016
{
    if (!is.null(design)) { 
        if (!missing(ncells)) { 
            stop("cannot specify both 'ncells' and 'design'")
        }
        ncells <- nrow(design) - qr(design)$rank
    }
    out <- .Call(cxx_get_null_rho, as.integer(ncells), as.integer(iters))
    if (is.character(out)) { 
        stop(out)
    }
    out <- sort(out)
    return(out)  
}

setGeneric("correlatePairs", function(x, ...) standardGeneric("correlatePairs"))

setMethod("correlatePairs", "ANY", function(x, null.dist=NULL, design=NULL, BPPARAM=bpparam(), use.names=TRUE)
# This calculates a (modified) Spearman's rho for each pair of genes.
#
# written by Aaron Lun
# created 10 February 2016
# last modified 21 February 2016
{
    exprs <- as.matrix(x)
    if (!is.null(design)) { 
        fit <- lm.fit(y=t(exprs), x=design)
        exprs <- t(fit$effects[-fit$qr$pivot[seq_len(fit$rank)],])
    }
    ncells <- ncol(exprs)
    if (is.null(null.dist)) { 
        null.dist <- correlateNull(ncells)
    } else {
        null.dist <- as.double(null.dist)
    }
    if (is.unsorted(null.dist)) { 
        null.dist <- sort(null.dist)
    }
    ranked.exprs <- apply(exprs, 1, FUN=rank, ties.method="random")

    # Generating all pairs of genes
    ngenes <- nrow(exprs)
    if (ngenes < 2L) { stop("need at least two genes to compute correlations") }
    all.pairs <- combn(ngenes, 2L)
    gene1 <- all.pairs[1,]
    gene2 <- all.pairs[2,]

    # Running through each set of jobs 
    workass <- .workerAssign(length(gene1), BPPARAM)
    out <- bplapply(seq_along(workass$start), FUN=function(core) {
        to.use <- workass$start[core]:workass$end[core]
        .Call(cxx_compute_rho, gene1[to.use], gene2[to.use], ncells, ranked.exprs, null.dist)
    }, BPPARAM=BPPARAM)

    # Peeling apart the output
    all.rho <- all.pval <- list()
    for (x in seq_along(out)) {
        current <- out[[x]]
        if (is.character(current)) { stop(current) }
        all.rho[[x]] <- current[[1]]
        all.pval[[x]] <- current[[2]]
    }
    all.pval <- unlist(all.pval)
    all.rho <- unlist(all.rho)

    # Returning some useful output
    newnames <- NULL
    if (is.logical(use.names)) {
        if (use.names) {
            newnames <- rownames(exprs)
        }
    } else if (is.character(use.names)) {
        if (length(use.names)!=nrow(exprs)) {
            stop("length of 'use.names' does not match 'exprs' nrow")
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

setMethod("correlatePairs", "SCESet", function(x, ..., assay="exprs", get.spikes=FALSE) {
    correlatePairs(.getUsedMatrix(x, assay, get.spikes), ...)             
})

