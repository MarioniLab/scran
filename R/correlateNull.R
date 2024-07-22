#' Build null correlations
#' 
#' Build a distribution of correlations under the null hypothesis of independent expression between pairs of genes.
#' This is now deprecated as \code{\link{correlatePairs}} uses an approximation instead.
#'
#' @param ncells An integer scalar indicating the number of cells in the data set.
#' @param iters An integer scalar specifying the number of values in the null distribution.
#' @param block A factor specifying the blocking level for each cell.
#' @param design A numeric design matrix containing uninteresting factors to be ignored.
#' @param equiweight A logical scalar indicating whether statistics from each block should be given equal weight.
#' Otherwise, each block is weighted according to its number of cells.
#' Only used if \code{block} is specified.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object that specifies the manner of parallel processing to use.
#'
#' @details
#' The \code{correlateNull} function constructs an empirical null distribution for Spearman's rank correlation when it is computed with \code{ncells} cells.
#' This is done by shuffling the ranks, calculating the correlation and repeating until \code{iters} values are obtained.
#' No consideration is given to tied ranks, which has implications for the accuracy of p-values in \code{\link{correlatePairs}}.
#'
#' If \code{block} is specified, a null correlation is created within each level of \code{block} using the shuffled ranks.
#' The final correlation is then defined as the average of the per-level correlations, 
#' weighted by the number of cells in that level if \code{equiweight=FALSE}.
#' Levels with fewer than 3 cells are ignored, and if no level has 3 or more cells, all returned correlations will be \code{NA}.
#'
#' If \code{design} is specified, the same process is performed on ranks derived from simulated residuals computed by fitting the linear model to a vector of normally distributed values.
#' If there are not at least 3 residual d.f., all returned correlations will be \code{NA}.
#' The \code{design} argument cannot be used at the same time as \code{block}.
#' 
#' % Yeah, we could use a t-distribution for this, but the empirical distribution is probably more robust if you have few cells (or effects, after batch correction).
#'
#' @return
#' A numeric vector of length \code{iters} is returned containing the sorted correlations under the null hypothesis of no correlations.
#'
#' @author
#' Aaron Lun
#' 
#' @seealso
#' \code{\link{correlatePairs}}, where the null distribution is used to compute p-values.
#'
#' @examples
#' set.seed(0)
#' ncells <- 100
#'
#' # Simplest case:
#' null.dist <- correlateNull(ncells, iters=10000)
#' hist(null.dist)
#'
#' # With a blocking factor:
#' block <- sample(LETTERS[1:3], ncells, replace=TRUE)
#' null.dist <- correlateNull(block=block, iters=10000)
#' hist(null.dist)
#'
#' # With a design matrix.
#' cov <- runif(ncells)
#' X <- model.matrix(~cov)
#' null.dist <- correlateNull(design=X, iters=10000)
#' hist(null.dist)
#' 
#' @export
#' @importFrom BiocParallel SerialParam bpmapply bpstart bpstop bpisup
#' @importFrom scuttle .bpNotSharedOrUp .ranksafeQR
correlateNull <- function(ncells, iters=1e6, block=NULL, design=NULL, equiweight=TRUE, BPPARAM=SerialParam()) {
    .Deprecated()
    if (!bpisup(BPPARAM)) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM))
    }

    iters <- as.integer(iters)
    iters.per.core <- .niters_by_nworkers(iters, BPPARAM)

    blockFUN <- function(group) {
        pcg.state <- .setup_pcg_state(iters.per.core)
        out.g <- bpmapply(Niters=iters.per.core, Seeds=pcg.state$seeds, Streams=pcg.state$streams,
            MoreArgs=list(Ncells=length(group)), FUN=get_null_rho,
            SIMPLIFY=FALSE, USE.NAMES=FALSE, BPPARAM=BPPARAM)
        unlist(out.g) 
    }

    designFUN <- function(design) {
        QR <- .ranksafeQR(design)
        pcg.state <- .setup_pcg_state(iters.per.core)
        out.g <- bpmapply(Niters=iters.per.core, Seeds=pcg.state$seeds, Streams=pcg.state$streams, 
            MoreArgs=list(qr=QR$qr, qraux=QR$qraux), FUN=get_null_rho_design,
            SIMPLIFY=FALSE, USE.NAMES=FALSE, BPPARAM=BPPARAM)
        unlist(out.g) 
    }

    output <- .correlator_base(ncells, block, design, equiweight, blockFUN, designFUN)
    if (is.null(output)) {
        output <- rep(NA_real_, iters)
    }
    output
}

.correlator_base <- function(ncells, block, design, equiweight, blockFUN, designFUN) {
    if (is.null(design)) { 
        if (is.null(block)) {
            block <- rep(1L, ncells)
        } else if (!missing(ncells) && ncells!=length(block)) { 
            stop("cannot specify both 'ncells' and 'block'")
        }

        groupings <- split(seq_along(block), block)
        if (equiweight) {
            weights <- rep(1, length(groupings))
        } else {
            weights <- lengths(groupings)
        }

        # Estimating the correlation as a weighted mean of the correlations in each group.
        # This avoids the need for the normality assumption in the residual effect simulation.
        out <- total <- 0
        for (g in seq_along(groupings)) {
            ngr <- length(groupings[[g]])
            if (ngr <= 2L) {
                next            
            }

            out.g <- blockFUN(groupings[[g]])
            w <- weights[g]
            out <- out + unlist(out.g) * w
            total <- total + w
        }
        out/total

    } else {
        if (!is.null(block)) {
            stop("cannot specify both 'block' and 'design'")
        }
        if (!missing(ncells) && ncells!=nrow(design)) {
            stop("cannot specify both 'ncells' and 'design'")
        }
        if (nrow(design) - ncol(design) > 2L) {
            designFUN(design)
        } else {
            NULL
        }
    }
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
