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
#' The correlation between a pair of genes is then computed from the residuals of the fitted model (see \code{\link{ResidualMatrix}} for the implementation details).
#' However, we recommend using \code{block} wherever possible as \code{design} assumes normality of errors and deals poorly with ties.
#' Specifically, zero counts within or across groups may no longer be tied when converted to residuals, potentially resulting in spuriously large correlations.
#'
#' If any level of \code{block} has fewer than 3 cells, it is ignored.
#' If all levels of \code{block} have fewer than 3 cells, all output statistics are set to \code{NA}.
#' Similarly, if \code{design} has fewer than 3 residual d.f., all output statistics are set to \code{NA}.
#'
#' @section Gene selection:
#' The \code{pairings} argument specifies the pairs of genes to compute correlations for:
#' \itemize{
#' \item By default, correlations will be computed between all pairs of genes with \code{pairings=NULL}.
#' Genes that occur earlier in \code{x} are labelled as \code{gene1} in the output DataFrame.
#' Redundant permutations are not reported.
#' \item If \code{pairings} is a list of two vectors, correlations will be computed between one gene in the first vector and another gene in the second vector.
#' This improves efficiency if the only correlations of interest are those between two pre-defined sets of genes.
#' Genes in the first vector are always reported as \code{gene1}.
#' \item If \code{pairings} is an integer/character matrix of two columns, each row is assumed to specify a gene pair.
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
#' @importFrom ResidualMatrix ResidualMatrix
#' @importFrom Matrix t
#' @importFrom BiocGenerics cbind
#' @importFrom metapod parallelStouffer
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

    # Checking which pairwise correlations should be computed.
    pair.out <- .construct_pair_indices(subset.row=subset.row, x=x, pairings=pairings)
    gene1 <- pair.out$gene1
    gene2 <- pair.out$gene2
    reorder <- pair.out$reorder

    if (!is.null(design)) {
        block <- NULL 
        x <- t(ResidualMatrix(t(x), design=design))
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
        output <- unlist(output)
        rhos[[b]] <- output

        p.out <- rhoToPValue(output, n=ncol(x0))
        up.p[[b]] <- p.out$positive
        down.p[[b]] <- p.out$negative
    }

    weights <- NULL
    if (!equiweight && length(block) > 1) {
        weights <- lengths(block)
    }
    all.rho <- .weighted_average_vals(rhos, weights)

    up.p <- parallelStouffer(up.p, weights=weights)$p.value
    down.p <- parallelStouffer(down.p, weights=weights)$p.value
    all.pval <- pmin(up.p, down.p) * 2

    if (!is.null(rownames(x)) && use.names) {
        gene1 <- rownames(x)[gene1]
        gene2 <- rownames(x)[gene2]
    } 

    out <- DataFrame(gene1=gene1, gene2=gene2, rho=all.rho, 
        p.value=all.pval, FDR=p.adjust(all.pval, method="BH"))
    if (reorder) {
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

#' @importFrom utils combn
#' @importFrom scuttle .subset2index
.construct_pair_indices <- function(subset.row, x, pairings) {
    subset.row <- .subset2index(subset.row, x, byrow=TRUE)
    reorder <- TRUE

    if (is.matrix(pairings)) {
        # If matrix, we're using pre-specified pairs.
        if ((!is.numeric(pairings) && !is.character(pairings)) || ncol(pairings)!=2L) { 
            stop("'pairings' should be a numeric/character matrix with 2 columns") 
        }
        s1 <- .subset2index(pairings[,1], x, byrow=TRUE)
        s2 <- .subset2index(pairings[,2], x, byrow=TRUE)

        # Discarding elements not in subset.row.
        keep <- s1 %in% subset.row & s2 %in% subset.row
        gene1 <- s1[keep]
        gene2 <- s2[keep]
        reorder <- FALSE

    } else if (is.list(pairings)) {
        # If list, we're correlating between one gene selected from each of two pools.
        if (length(pairings)!=2L) { 
            stop("'pairings' as a list should have length 2") 
        }
        converted <- lapply(pairings, FUN=function(gene.set) {
            gene.set <- .subset2index(gene.set, x, byrow=TRUE)
            intersect(gene.set, subset.row) # automatically gets rid of duplicates.
        })
        if (any(lengths(converted)==0L)) { 
            stop("need at least one gene in each set to compute correlations") 
        }

        all.pairs <- expand.grid(converted[[1]], converted[[2]])
        keep <- all.pairs[,1]!=all.pairs[,2]
        gene1 <- all.pairs[keep,1]
        gene2 <- all.pairs[keep,2]

    } else if (is.null(pairings)) {
        # Otherwise, it's assumed to be a single pool, and we're just correlating between pairs within it.
        ngenes <- length(subset.row)
        if (ngenes < 2L) { 
            stop("need at least two genes to compute correlations") 
        }
       
        # Generating all pairs of genes within the subset.
        all.pairs <- combn(subset.row, 2L)
        gene1 <- all.pairs[1,]
        gene2 <- all.pairs[2,]

    } else {
        stop("pairings should be a list, matrix or NULL")
    }

    list(gene1=gene1, gene2=gene2, reorder=reorder)
}

.choose_gene_names <- function(x, use.names) {
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
    if (any(results$FDR > threshold & results$limited, na.rm=TRUE)) { 
        warning(sprintf("lower bound on p-values at a FDR of %s, increase 'iter'", as.character(threshold)))
    }
    invisible(NULL)
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
