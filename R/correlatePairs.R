#' Test for significant correlations
#' 
#' Identify pairs of genes that are significantly correlated in their expression profiles, based on Spearman's rank correlation.
#'
#' @param design A numeric design matrix containing uninteresting factors to be ignored.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object that specifies the manner of parallel processing to use.
#' @param x 
#'     A numeric matrix-like object of log-normalized expression values, where rows are genes and columns are cells.
#'     Alternatively, a \linkS4class{SummarizedExperiment} object containing such a matrix.
#' 
#' @param null.dist,ties.method,iters Deprecated arguments, ignored.
#' @param use.names 
#'     A logical scalar specifying whether the row names of \code{x} should be used in the output.
#'     Alternatively, a character vector containing the names to use.
#' @param subset.row See \code{?"\link{scran-gene-selection}"}.
#' @param pairings A \code{NULL} value indicating that all pairwise correlations should be computed;
#' or a list of 2 vectors of genes between which correlations are to be computed;
#' or a integer/character matrix with 2 columns of specific gene pairs - see below for details.
#' @param block A factor specifying the blocking level for each cell in \code{x}.
#' If specified, correlations are computed separately within each block and statistics are combined across blocks.
#' @param equiweight A logical scalar indicating whether statistics from each block should be given equal weight.
#' Otherwise, each block is weighted according to its number of cells.
#' Only used if \code{block} is specified.
#' @param ... 
#' For the generic, additional arguments to pass to specific methods.
#'
#' For the SummarizedExperiment method, additional methods to pass to the ANY method.
#' @param assay.type A string specifying which assay values to use.
#' 
#' @details
#' The \code{correlatePairs} function identifies significant correlations between all pairs of genes in \code{x}.
#' This allows prioritization of genes that are driving systematic substructure in the data set.
#' By definition, such genes should be correlated as they are behaving in the same manner across cells.
#' In contrast, genes driven by random noise should not exhibit any correlations with other genes.
#' 
#' We use Spearman's rho to quantify correlations robustly based on ranked gene expression.
#' To identify correlated gene pairs, the significance of non-zero correlations is assessed using \code{\link{rhoToPValue}}.
#' The null hypothesis is that the ranking of normalized expression across cells should be independent between genes.
#' Correction for multiple testing is done using the BH method.
#'
#' For the SingleCellExperiment method, normalized expression values should be specified by \code{assay.type}. 
#' 
#' @return
#' A \linkS4class{DataFrame} is returned with one row per gene pair and the following fields:
#' \describe{
#' \item{\code{gene1, gene2}:}{
#'     Character or integer fields specifying the genes in the pair.
#'     If \code{use.names=FALSE}, integers are returned representing row indices of \code{x}, otherwise gene names are returned.
#' }
#' \item{\code{rho}:}{A numeric field containing the approximate Spearman's rho.}
#' \item{\code{p.value, FDR}:}{Numeric fields containing the approximate p-value and its BH-corrected equivalent.}
#' } 
#' Rows are sorted by increasing \code{p.value} and, if tied, decreasing absolute size of \code{rho}.
#' The exception is if \code{subset.row} is a matrix, in which case each row in the dataframe correspond to a row of \code{subset.row}.
#' 
#' @section Accounting for uninteresting variation:
#' If the experiment has known (and uninteresting) factors of variation, these can be included in \code{design} or \code{block}.
#' \code{correlatePairs} will then attempt to ensure that these factors do not drive strong correlations between genes.
#' Examples might be to block on batch effects or cell cycle phase, which may have substantial but uninteresting effects on expression.
#' 
#' The approach used to remove these factors depends on whether \code{design} or \code{block} is used.
#' If there is only one factor, e.g., for plate or animal of origin, \code{block} should be used.
#' Each level of the factor is defined as a separate group of cells.
#' For each pair of genes, correlations are computed within each group, and a mean of the correlations is taken across all groups.
#' If \code{equiweight=FALSE}, a weighted mean is computed based on the size of each group.
#'
#' Similarly, \code{\link{parallelStouffer}} is used to combine the (one-sided) p-values across all groups.
#' This is done for each direction and a final p-value is computed for each gene pair using this Bonferri method.
#' The idea is to ensure that the final p-value is only low when correlations are in the same direction across groups.
#' If \code{equiweight=FALSE}, each p-value is weighted by the size of the corresponding group.
#' 
#' For experiments containing multiple factors or covariates, a design matrix should be passed into \code{design}.
#' The correlation between each pair of genes is then computed from the residuals of the fitted model.
#' However, we recommend using \code{block} wherever possible as \code{design} assumes normality of errors and deals poorly with ties.
#' Specifically, zero counts within or across groups may no longer be tied when converted to residuals, potentially resulting in spuriously large correlations.
#'
#' If any level of \code{block} has fewer than 3 cells, it is ignored.
#' If all levels of \code{block} have fewer than 3 cells, all output statistics are set to \code{NA}.
#' Similarly, if \code{design} has fewer than 3 residual d.f., all output statistics are set to \code{NA}.
#'
#' @section Gene selection:
#' The \code{pairings} argument specifies the pairs of genes that should be used to compute correlations.
#' This can be:
#' \itemize{
#' \item \code{NULL}, in which case correlations will be computed between all pairs of genes in \code{x}.
#' Genes that occur earlier in \code{x} are labelled as \code{gene1} in the output DataFrame.
#' Redundant permutations are not reported.
#' \item A list of two vectors, where each list element defines a subset of genes in \code{x} as an integer, character or logical vector.
#' In this case, correlations will be computed between one gene in the first vector and another gene in the second vector.
#' This improves efficiency if the only correlations of interest are those between two pre-defined sets of genes.
#' Genes in the first vector are always reported as \code{gene1}.
#' \item An integer/character matrix of two columns.
#' In this case, each row is assumed to specify a gene pair based on the row indices (integer) or row names (character) of \code{x}.
#' Correlations will then be computed for only those gene pairs, and the returned dataframe will \emph{not} be sorted by p-value.
#' Genes in the first column of the matrix are always reported as \code{gene1}.
#' }
#' 
#' If \code{subset.row} is not \code{NULL}, only the genes in the selected subset are used to compute correlations - see \code{?"\link{scran-gene-selection}"}.
#' This will interact properly with \code{pairings}, such that genes in \code{pairings} and not in \code{subset.row} will be ignored.
#' 
#' We recommend setting  \code{subset.row} and/or \code{pairings} to contain only the subset of genes of interest.
#' This reduces computational time and memory usage by only computing statistics for the gene pairs of interest.
#' For example, we could select only HVGs to focus on genes contributing to cell-to-cell heterogeneity (and thus more likely to be involved in driving substructure).
#' There is no need to account for HVG pre-selection in multiple testing, because rank correlations are unaffected by the variance.
#' 
#' @author
#' Aaron Lun
#' 
#' @seealso
#' Compare to \code{\link{cor}} for the standard Spearman's calculation.
#' 
#' Use \code{\link{correlateGenes}} to get per-gene correlation statistics.
#' 
#' @references
#' Lun ATL (2019).
#' Some thoughts on testing for correlations.
#' \url{https://ltla.github.io/SingleCellThoughts/software/correlations/corsim.html}
#' 
#' @examples
#' library(scuttle)
#' sce <- mockSCE()
#' sce <- logNormCounts(sce)
#'
#' # Basic pairwise application.
#' out <- correlatePairs(sce, subset.row=1:100)
#' head(out)
#'
#' # Computing between specific subsets of genes:
#' out <- correlatePairs(sce, pairings=list(1:10, 110:120))
#' head(out)
#'
#' # Computing between specific pairs:
#' out <- correlatePairs(sce, pairings=rbind(c(1,10), c(2, 50)))
#' head(out)
#' 
#' @name correlatePairs
NULL

#' @importFrom BiocParallel SerialParam bpstart bpstop
#' @importFrom S4Vectors DataFrame
#' @importFrom stats p.adjust
#' @importFrom scuttle .bpNotSharedOrUp .assignIndicesToWorkers .splitVectorByWorkers fitLinearModel
#' @importFrom DelayedArray DelayedArray
#' @importFrom Matrix t
#' @importFrom BiocGenerics cbind
#' @importFrom metapod parallelStouffer averageParallelStats
.correlate_pairs <- function(x, null.dist=NULL, ties.method=NULL, iters=NULL, 
    block=NULL, design=NULL, equiweight=TRUE, use.names=TRUE, subset.row=NULL, 
    pairings=NULL, BPPARAM=SerialParam())
{
    if (.bpNotSharedOrUp(BPPARAM)) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM))
    }

    if (!is.null(iters) || !is.null(ties.method) || !is.null(null.dist)) {
        .Deprecated(msg="'iters=', 'ties.method=' and 'null.dist=' are deprecated.")
    }

    pair.out <- .process_corr_pairings(pairings, x, subset.row=subset.row)
    gene1 <- pair.out$id1
    gene2 <- pair.out$id2

    if (!is.null(design)) {
        block <- NULL 
        x <- t(ResidualMatrix::ResidualMatrix(t(x), design=design))
    }

    if (is.null(block)) {
        block <- list(NULL)
    } else {
        block <- split(seq_along(block), block)
    }
    rhos <- up.p <- down.p <- vector("list", length(block)) 

    x <- DelayedArray(x)
    for (b in seq_along(block)) {
        x0 <- x
        if (!is.null(block[[b]])) {
            x0 <- x0[,block[[b]],drop=FALSE]
        }

        x1 <- x0[gene1,,drop=FALSE]
        x2 <- x0[gene2,,drop=FALSE]
        y <- cbind(x1, x2) # so we can block process both matrices at the same time.

        output <- rowBlockApply(y, FUN=.calculate_rho, BPPARAM=BPPARAM)
        output <- unlist(output, use.names=FALSE)
        rhos[[b]] <- unname(output)

        p.out <- rhoToPValue(output, n=ncol(x0))
        up.p[[b]] <- p.out$positive
        down.p[[b]] <- p.out$negative
    }

    weights <- NULL
    if (!equiweight && length(block) > 1) {
        weights <- lengths(block)
    }
    all.rho <- averageParallelStats(rhos, weights)

    up.p <- parallelStouffer(up.p, weights=weights)$p.value
    down.p <- parallelStouffer(down.p, weights=weights)$p.value
    all.pval <- pmin(up.p, down.p) * 2

    if (!is.null(rownames(x)) && use.names) {
        gene1 <- rownames(x)[gene1]
        gene2 <- rownames(x)[gene2]
    } 

    out <- DataFrame(gene1=gene1, gene2=gene2, rho=all.rho, p.value=all.pval, FDR=p.adjust(all.pval, method="BH"))

    if (pair.out$mode != "predefined pairs") {
        # Preserving the order if the pairings were predefined as a matrix.
        out <- out[order(out$p.value, -abs(out$rho)),]
        rownames(out) <- NULL
    }

    out
}

###############################
### INTERNAL (calculations) ###
###############################

#' @importFrom BiocParallel bpmapply 
#' @importFrom DelayedMatrixStats rowRanks rowSds
#' @importFrom Matrix t rowMeans
#' @importFrom stats var
.calculate_rho <- function(block) {
    half <- ncol(block)/2
    if (half < 2) {
        return(rep(NA_real_, nrow(block)))
    }

    first <- seq_len(half)
    halves <- list(block[,first,drop=FALSE], block[,half + first,drop=FALSE])

    for (i in seq_along(halves)) {
        current <- halves[[i]]
        current <- rowRanks(current, ties.method="average")
        current <- (current - rowMeans(current)) / rowSds(current)
        halves[[i]] <- current
    }

    out <- rowSums(halves[[1]] * halves[[2]]) / (half - 1)
    out[!is.finite(out)] <- NA_real_
    out
}

#' @importFrom scuttle .subset2index
.process_corr_pairings <- function(pairings, x, subset.row=NULL) {
    if (!is.null(subset.row)) {
        subset.row <- .subset2index(subset.row, x, byrow=TRUE)
    }

    # Unlike .subset2index, this doesn't throw on NAs.
    dummy <- seq_len(nrow(x))
    names(dummy) <- rownames(x)

    .SUBSET <- function(request, clean=TRUE) {
        if (is.null(request)) {
            if (!is.null(subset.row)) {
                out <- subset.row
            } else {
                out <- unname(dummy)
            }
        } else {
            out <- unname(dummy[request])
            if (!is.null(subset.row)) {
                out[!out %in% subset.row] <- NA_integer_
            }
        }
        if (clean) {
            out <- unique(out[!is.na(out)])
        }
        out
    }

    pair.out <- .expand_pairings_core(pairings, .SUBSET)

    if (pair.out$mode == "single pool") {
        # Removing redundant permutations, if we were asked to generate them from a single pool.
        keep <- pair.out$id1 < pair.out$id2
        pair.out$id1 <- pair.out$id1[keep]
        pair.out$id2 <- pair.out$id2[keep]
    }

    pair.out
}

#############################
### INTERNAL (S4 methods) ###
#############################

#' @export
#' @rdname correlatePairs
setGeneric("correlatePairs", function(x, ...) standardGeneric("correlatePairs"))

#' @export
#' @rdname correlatePairs
setMethod("correlatePairs", "ANY", .correlate_pairs)

#' @export
#' @rdname correlatePairs
#' @importFrom SummarizedExperiment assay
setMethod("correlatePairs", "SummarizedExperiment", function(x, ..., assay.type="logcounts") {
    .correlate_pairs(assay(x, i=assay.type), ...)
})
