#' @export
#' @importFrom BiocParallel SerialParam bpmapply
correlateNull <- function(ncells, iters=1e6, block=NULL, design=NULL, BPPARAM=SerialParam()) 
# This builds a null distribution for the modified Spearman's rho.
#
# written by Aaron Lun
# created 10 February 2016    
{
    if (!is.null(block)) { 
        if (!missing(ncells)) { 
            stop("cannot specify both 'ncells' and 'block'")
        }
        groupings <- table(block)

        # Estimating the correlation as a weighted mean of the correlations in each group.
        # This avoids the need for the normality assumption in the residual effect simulation.
        out <- 0
        for (ngr in groupings) {
            out.g <- .within_block_null(ncells=ngr, iters=as.integer(iters), BPPARAM=BPPARAM)
            out <- out + out.g * ngr
        }
        out <- out/length(block)
        attrib <- list(block=block)
        
    } else if (!is.null(design)) { 
        if (!missing(ncells)) { 
            stop("cannot specify both 'ncells' and 'design'")
        }

        # Using residual effects to compute the correlations.
        QR <- .ranksafe_qr(design)
        iters.per.core <- .niters_by_nworkers(as.integer(iters), BPPARAM)
        seeds.per.core <- lapply(iters.per.core, .create_seeds)
        out <- bpmapply(iters=iters.per.core, seeds=seeds.per.core, MoreArgs=list(QR=QR), FUN=.with_design_null, SIMPLIFY=FALSE, USE.NAMES=FALSE)
        out <- unlist(out)
        attrib <- list(design=design)

    } else {
        out <- .within_block_null(iters=as.integer(iters), ncells=as.integer(ncells), BPPARAM=BPPARAM)
        attrib <- NULL
    }

    # Storing attributes, to make sure it matches up.
    out <- sort(out)
    attributes(out) <- attrib
    return(out)  
}

#### Internal functions (no design) ####

#' @importFrom BiocParallel bpmapply
.within_block_null <- function(iters, ncells, BPPARAM) {
    iters.per.core <- .niters_by_nworkers(as.integer(iters), BPPARAM)
    seeds.per.core <- lapply(iters.per.core, .create_seeds)
    out <- bpmapply(iters=iters.per.core, seeds=seeds.per.core, MoreArgs=list(ncells=as.integer(ncells)), 
        FUN=.no_design_null, SIMPLIFY=FALSE, USE.NAMES=FALSE)
    unlist(out)
}

.no_design_null <- function(ncells, iters, seeds) {
    .Call(cxx_get_null_rho, ncells, iters, seeds)
}

.with_design_null <- function(iters, QR, seeds) {
    .Call(cxx_get_null_rho_design, QR$qr, QR$qraux, iters, seeds)
}

#' @importFrom BiocParallel bpnworkers
.niters_by_nworkers <- function(iters, BPPARAM) {
    nworkers <- bpnworkers(BPPARAM)
    if (iters <= nworkers) {
        jobs <- integer(nworkers)
        jobs[seq_len(iters)] <- 1L
    } else {
        per.worker <- as.integer(floor(iters/nworkers))
        jobs <- rep(per.worker, nworkers)
        jobs[1] <- iters - sum(jobs[-1])
    }
    jobs
}
