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
#' @param null.dist A numeric vector of rho values under the null hypothesis.
#' @param ties.method String specifying how tied ranks should be handled.
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
#' @param iters Integer scalar specifying the number of iterations to use in \code{\link{correlateNull}} to build the null distribution.
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
#' To identify correlated gene pairs, the significance of non-zero correlations is assessed using a permutation test.
#' The null hypothesis is that the ranking of normalized expression across cells should be independent between genes.
#' This allows us to construct a null distribution by randomizing the ranks within each gene.
#' A pre-computed empirical distribution can be supplied as \code{null.dist}, which saves some time by avoiding repeated calls to \code{\link{correlateNull}}.
#' 
#' The p-value for each gene pair is defined as the tail probability of this distribution at the observed correlation (with some adjustment to avoid zero p-values).
#' Correction for multiple testing is done using the BH method.
#' The lower bound of the p-values is determined by the value of \code{iters}.
#' If the \code{limited} field is \code{TRUE} in the returned dataframe, it may be possible to obtain lower p-values by increasing \code{iters}.
#' This should be examined for non-significant pairs, in case some correlations are overlooked due to computational limitations.
#' The function will automatically raise a warning if any genes are limited in their significance at a FDR of 5\%.
#'
#' For the SingleCellExperiment method, correlations should be computed for normalized expression values in the specified \code{assay.type}. 
#' 
#' @return
#' A \linkS4class{DataFrame} is returned with one row per gene pair and the following fields:
#' \describe{
#' \item{\code{gene1, gene2}:}{
#'     Character or integer fields specifying the genes in the pair.
#'     If \code{use.names=FALSE}, integers are returned representing row indices of \code{x}, otherwise gene names are returned.
#' }
#' \item{\code{rho}:}{A numeric field containing the approximate Spearman's rho.}
#' \item{\code{p.value, FDR}:}{Numeric fields containing the permutation p-value and its BH-corrected equivalent.}
#' \item{\code{limited}:}{A logical scalar indicating whether the p-value is at its lower bound, defined by \code{iters}.}
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
#' For each pair of genes, correlations are computed within each group, and a weighted mean based on the group size) of the correlations is taken across all groups.
#' The same strategy is used to generate the null distribution where ranks are computed and shuffled within each group.
#' 
#' For experiments containing multiple factors or covariates, a design matrix should be passed into \code{design}.
#' The \code{correlatePairs} function will fit a linear model to the (log-normalized) expression values.
#' The correlation between a pair of genes is then computed from the residuals of the fitted model.
#' Similarly, to obtain a null distribution of rho values, normally-distributed random errors are simulated in a fitted model based on \code{design};
#'     the corresponding residuals are generated from these errors; and the correlation between sets of residuals is computed at each iteration.
#' 
#' We recommend using \code{block} wherever possible.
#' While \code{design} can also be used for one-way layouts, this is not ideal as it assumes normality of errors and deals poorly with ties.
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
#' @section Handling tied values:
#' By default, \code{ties.method="expected"} which uses the expectation of the pairwise correlation from randomly broken ties.
#' Imagine obtaining unique ranks by randomly breaking ties in expression values, e.g., by addition of some very small random jitter.
#' The reported value of the correlation is simply the expected value across all possible permutations of broken ties.
#' 
#' With \code{ties.method="average"}, each set of tied observations is assigned the average rank across all tied values.
#' This yields the standard value of Spearman's rho as computed by \code{\link{cor}}.
#' 
#' We use the expected rho by default as avoids inflated correlations in the presence of many ties.
#' This is especially relevant for single-cell data containing many zeroes that remain tied after scaling normalization.
#' An extreme example is that of two genes that are only expressed in the same cell, for which the standard rho is 1 while the expected rho is close to zero.
#' 
#' Note that the p-value calculations are not accurate in the presence of ties, as tied ranks are not considered by \code{correlateNull}.
#' With \code{ties.method="expected"}, the p-values are probably a bit conservative.
#' With \code{ties.method="average"}, they will be very anticonservative.
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
#' Phipson B and Smyth GK (2010).
#' Permutation P-values should never be zero: calculating exact P-values when permutations are randomly drawn.
#' \emph{Stat. Appl. Genet. Mol. Biol.} 9:Article 39.
#' 
#' @examples
#' library(scater)
#' sce <- mockSCE()
#' sce <- logNormCounts(sce)
#'
#' # Basic pairwise application (turning down iters for speed).
#' out <- correlatePairs(sce, subset.row=1:100, iters=1e5)
#' head(out)
#'
#' # Computing between specific subsets of genes:
#' out <- correlatePairs(sce, pairings=list(1:10, 110:120), iters=1e5)
#' head(out)
#'
#' # Computing between specific pairs:
#' out <- correlatePairs(sce, pairings=rbind(c(1,10), c(2, 50)), iters=1e5)
#' head(out)
#' 
#' @name correlatePairs
NULL

#' @importFrom BiocParallel SerialParam bpstart bpstop
#' @importFrom S4Vectors DataFrame
#' @importFrom stats p.adjust
#' @importFrom scuttle .bpNotSharedOrUp .assignIndicesToWorkers .splitVectorByWorkers fitLinearModel
.correlate_pairs <- function(x, null.dist=NULL, ties.method=c("expected", "average"), 
    iters=1e6, block=NULL, design=NULL, equiweight=TRUE, use.names=TRUE, subset.row=NULL, 
    pairings=NULL, BPPARAM=SerialParam())
{
    # Checking which pairwise correlations should be computed.
    pair.out <- .construct_pair_indices(subset.row=subset.row, x=x, pairings=pairings)
    subset.row <- pair.out$subset.row
    gene1 <- pair.out$gene1
    gene2 <- pair.out$gene2
    reorder <- pair.out$reorder

    if (.bpNotSharedOrUp(BPPARAM)) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM))
    }

    null.dist <- .check_null_dist(x, block=block, design=design, iters=iters, 
        equiweight=equiweight, null.dist=null.dist, BPPARAM=BPPARAM)
    ties.method <- match.arg(ties.method)

    # Splitting up gene pairs into jobs for multicore execution, converting to 0-based indices.
    wout <- .assignIndicesToWorkers(length(gene1), BPPARAM)
    sgene1 <- .splitVectorByWorkers(gene1 - 1L, BPPARAM, assignments=wout)
    sgene2 <- .splitVectorByWorkers(gene2 - 1L, BPPARAM, assignments=wout)

    blockFUN <- function(subset.col) {
        .calculate_rho(sgene1, sgene2, x=x, subset.row=subset.row, 
            subset.col=subset.col, ties.method=ties.method, BPPARAM=BPPARAM)
    }

    designFUN <- function(design) {
        fitted <- fitLinearModel(x, design, subset.row=subset.row)
        resids <- x[subset.row,,drop=FALSE] - tcrossprod(fitted$coefficients, design)
        .calculate_rho(sgene1, sgene2, x=resids, subset.row=NULL, 
            subset.col=NULL, ties.method=ties.method, BPPARAM=BPPARAM)
    }

    all.rho <- .correlator_base(ncol(x), block, design, equiweight, blockFUN, designFUN)
    if (is.null(all.rho)) {
        all.rho <- rep(NA_real_, length(gene1))
    }

    # Computing p-values and formatting the output.
    stats <- .rho_to_pval(all.rho, null.dist)
    all.pval <- stats$p
    all.lim <- stats$limited

    final.names <- .choose_gene_names(subset.row=subset.row, x=x, use.names=use.names)
    gene1 <- final.names[gene1]
    gene2 <- final.names[gene2]

    out <- DataFrame(gene1=gene1, gene2=gene2, rho=all.rho, p.value=all.pval, 
        FDR=p.adjust(all.pval, method="BH"), limited=all.lim)
    if (reorder) {
        out <- out[order(out$p.value, -abs(out$rho)),]
        rownames(out) <- NULL
    }
    .is_sig_limited(out)

    out
}

##########################################
### INTERNAL (correlation calculation) ###
##########################################

.check_null_dist <- function(x, block, design, ..., null.dist) 
# This makes sure that the null distribution is in order.
{
    if (is.null(null.dist)) { 
        if (!is.null(block)) { 
            null.dist <- correlateNull(block=block, ...)
        } else if (!is.null(design)) { 
            null.dist <- correlateNull(design=design, ...)
        } else {
            null.dist <- correlateNull(ncol(x), ...)
        }
    }

    null.dist <- as.double(null.dist)
    if (is.unsorted(null.dist)) { 
        null.dist <- sort(null.dist)
    }
    null.dist
}

#' @importFrom BiocParallel bpmapply 
#' @importFrom DelayedMatrixStats rowRanks rowVars
#' @importFrom DelayedArray DelayedArray rowMeans getAutoBPPARAM setAutoBPPARAM
#' @importFrom Matrix t
#' @importFrom stats var
.calculate_rho <- function(sgene1, sgene2, x, subset.row, subset.col, ties.method, BPPARAM)
# Iterating through all blocking levels (for one-way layouts; otherwise, this is a loop of length 1).
# Computing correlations between gene pairs, and adding a weighted value to the final average.
{
    old <- getAutoBPPARAM()
    setAutoBPPARAM(BPPARAM)
    on.exit(setAutoBPPARAM(old))

    ranks <- rowRanks(DelayedArray(x), rows=subset.row, cols=subset.col, ties.method="average") 
    ranks <- DelayedArray(ranks)
    ranks <- ranks - rowMeans(ranks)

    if (ties.method=="average") {
        rank.scale <- rowVars(ranks)
    } else {
        rank.scale <- var(seq_len(ncol(ranks)))
    }

    N <- ncol(ranks)
    rank.scale <- rank.scale * (N-1)/N # convert from var to mean square.
    ranks <- ranks/sqrt(rank.scale)

    # Transposing for easier C++ per-gene access.
    # Realizing to avoid need to cache repeatedly.
    ranks <- t(ranks)
    ranks <- as.matrix(ranks) 

    out <- bpmapply(FUN=compute_rho_pairs, gene1=sgene1, gene2=sgene2, 
        MoreArgs=list(ranks=ranks), BPPARAM=BPPARAM, SIMPLIFY=FALSE)
    unlist(out)
}

.rho_to_pval <- function(all.rho, null.dist) 
# Estimating the p-values (need to shift values to break ties conservatively by increasing the p-value).
{
    left <- findInterval(all.rho + 1e-8, null.dist)
    right <- length(null.dist) - findInterval(all.rho - 1e-8, null.dist)
    limited <- left==0L | right==0L
    all.pval <- (pmin(left, right)+1)*2/(length(null.dist)+1)
    all.pval <- pmin(all.pval, 1)
    list(p=all.pval, limited=limited)
}

##################################
### INTERNAL (pair definition) ###
##################################

#' @importFrom utils combn
#' @importFrom scuttle .subset2index
.construct_pair_indices <- function(subset.row, x, pairings) 
# This returns a new subset-by-row vector, along with the pairs of elements
# indexed along that vector (i.e., "1" refers to the first element of subset.row,
# rather than the first element of "x").
{
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
        s1 <- s1[keep]
        s2 <- s2[keep]

        subset.row <- sort(unique(c(s1, s2)))
        gene1 <- match(s1, subset.row)
        gene2 <- match(s2, subset.row)
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

        subset.row <- sort(unique(unlist(converted)))
        m1 <- match(converted[[1]], subset.row)
        m2 <- match(converted[[2]], subset.row)
        all.pairs <- expand.grid(m1, m2)

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
        all.pairs <- combn(ngenes, 2L)
        gene1 <- all.pairs[1,]
        gene2 <- all.pairs[2,]

    } else {
        stop("pairings should be a list, matrix or NULL")
    }

    return(list(subset.row=subset.row, gene1=gene1, gene2=gene2, reorder=reorder))
}

####################################
### INTERNAL (output formatting) ###
####################################

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
