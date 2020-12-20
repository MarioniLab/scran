#' Cell cycle phase classification
#' 
#' Classify single cells into their cell cycle phases based on gene expression data.
#' 
#' @param x A numeric matrix-like object of gene expression values where rows are genes and columns are cells.
#' 
#' Alternatively, a \linkS4class{SummarizedExperiment} object containing such a matrix.
#' @param pairs A list of data.frames produced by \code{\link{sandbag}}, containing pairs of marker genes.
#' @param gene.names A character vector of gene names, with one value per row in \code{x}.
#' @param iter An integer scalar specifying the number of iterations for random sampling to obtain a cycle score.
#' @param min.iter An integer scalar specifying the minimum number of iterations for score estimation.
#' @param min.pairs An integer scalar specifying the minimum number of pairs for cycle estimation.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object to use for parallel processing across cells.
#' @param verbose A logical scalar specifying whether diagnostics should be printed to screen.
#' @param subset.row See \code{?"\link{scran-gene-selection}"}.
#' @param ... For the generic, additional arguments to pass to specific methods.
#'
#' For the SummarizedExperiment method, additional arguments to pass to the ANY method.
#' @param assay.type A string specifying which assay values to use, e.g., \code{"counts"} or \code{"logcounts"}.
#' 
#' @details
#' This function implements the classification step of the pair-based prediction method described by Scialdone et al. (2015).
#' To illustrate, consider classification of cells into G1 phase.
#' Pairs of marker genes are identified with \code{\link{sandbag}}, where the expression of the first gene in the training data is greater than the second in G1 phase but less than the second in all other phases.
#' For each cell, \code{cyclone} calculates the proportion of all marker pairs where the expression of the first gene is greater than the second in the new data \code{x} (pairs with the same expression are ignored).
#' A high proportion suggests that the cell is likely to belong in G1 phase, as the expression ranking in the new data is consistent with that in the training data.
#' 
#' Proportions are not directly comparable between phases due to the use of different sets of gene pairs for each phase.
#' Instead, proportions are converted into scores (see below) that account for the size and precision of the proportion estimate. 
#' The same process is repeated for all phases, using the corresponding set of marker pairs in \code{pairs}.
#' Cells with G1 or G2M scores above 0.5 are assigned to the G1 or G2M phases, respectively.
#' (If both are above 0.5, the higher score is used for assignment.)
#' Cells can be assigned to S phase based on the S score, but a more reliable approach is to define S phase cells as those with G1 and G2M scores below 0.5.
#'
#' Pre-trained classifiers are provided for mouse and human datasets, see \code{?\link{sandbag}} for more details.
#' However, note that the classifier may not be accurate for data that are substantially different from those used in the training set, e.g., due to the use of a different protocol.
#' In such cases, users can construct a custom classifier from their own training data using the \code{\link{sandbag}} function.
#' This is usually necessary for other model organisms where pre-trained classifiers are not available.
#' 
#' Users should \emph{not} filter out low-abundance genes before applying \code{cyclone}.
#' Even if a gene is not expressed in any cell, it may still be useful for classification if it is phase-specific.
#' Its lack of expression relative to other genes will still yield informative pairs, and filtering them out would reduce power.
#'
#' @section Description of the score calculation:
#' To make the proportions comparable between phases, a distribution of proportions is constructed by shuffling the expression values within each cell and recalculating the proportion.
#' The phase score is defined as the lower tail probability at the observed proportion.
#' High scores indicate that the proportion is greater than what is expected by chance if the expression of marker genes were independent 
#' (i.e., with no cycle-induced correlations between marker pairs within each cell).
#' 
#' % The shuffling assumes that the marker genes are IID from the same distribution of expression values, such that there's no correlations.
#' % The question is then what distribution of expression values to use - see below.
#' % Training also should protect against non-cycle-based correlations, as such they should be present across all phases and not get included in the marker set.
#' 
#' By default, shuffling is performed \code{iter} times to obtain the distribution from which the score is estimated.
#' However, some iterations may not be used if there are fewer than \code{min.pairs} pairs with different expression, such that the proportion cannot be calculated precisely.
#' A score is only returned if the distribution is large enough for stable calculation of the tail probability, i.e., consists of results from at least \code{min.iter} iterations.
#' 
#' Note that the score calculation in \code{cyclone} is slightly different from that described originally by Scialdone et al.
#' The original code shuffles all expression values within each cell, while in this implementation, only the expression values of genes in the marker pairs are shuffled.
#' This modification aims to use the most relevant expression values to build the null score distribution.
#' 
#' % In theory, this shouldn't matter, as the score calculation depends on the ranking of each gene.
#' % That should be the same regardless of the distribution of expression values -- each set of rankings is equally likely, no matter what.
#' % In practice, the number of tied expression values will differ between different set of genes, e.g., due to abundance (low counts more likely to get ties).
#' % The most appropriate comparison would involve the same number of ties as that used to calculate the observed score.
#' % It doesn't make sense, for example, to shuffle in a whole bunch of non-expressed genes (lots of zeroes, ties) when the markers are always expressed.
#' 
#' @return
#' A list is returned containing:
#' \describe{
#' \item{\code{phases}:}{A character vector containing the predicted phase for each cell.} 
#' \item{\code{scores}:}{A data frame containing the numeric phase scores for each phase and cell (i.e., each row is a cell).}
#' \item{\code{normalized.scores}:}{A data frame containing the row-normalized scores (i.e., where the row sum for each cell is equal to 1).}
#' }
#' 
#' @author
#' Antonio Scialdone,
#' with modifications by Aaron Lun
#' 
#' @seealso
#' \code{\link{sandbag}}, to generate the pairs from reference data.
#' 
#' @examples
#' set.seed(1000)
#' library(scuttle)
#' sce <- mockSCE(ncells=200, ngenes=1000)
#' 
#' # Constructing a classifier:
#' is.G1 <- which(sce$Cell_Cycle %in% c("G1", "G0"))
#' is.S <- which(sce$Cell_Cycle=="S")
#' is.G2M <- which(sce$Cell_Cycle=="G2M")
#' out <- sandbag(sce, list(G1=is.G1, S=is.S, G2M=is.G2M))
#' 
#' # Classifying a new dataset:
#' test <- mockSCE(ncells=50)
#' assignments <- cyclone(test, out)
#' head(assignments$scores)
#' table(assignments$phases)
#' 
#' @references
#' Scialdone A, Natarajana KN, Saraiva LR et al. (2015). 
#' Computational assignment of cell-cycle stage from single-cell transcriptome data.
#' \emph{Methods} 85:54--61
#'
#' @name cyclone
NULL

#' @export
#' @rdname cyclone
setGeneric("cyclone", function(x, ...) standardGeneric("cyclone"))

#' @importFrom BiocParallel SerialParam bpstart bpstop
#' @importFrom scuttle .bpNotSharedOrUp .subset2index 
#' @importFrom beachmat colBlockApply
.cyclone <- function(x, pairs, gene.names=rownames(x), iter=1000, min.iter=100, min.pairs=50, 
    BPPARAM=SerialParam(), verbose=FALSE, subset.row=NULL)
{ 
    if (length(gene.names)!=nrow(x)) {
        stop("length of 'gene.names' must be equal to 'x' nrows")
    }
    iter <- as.integer(iter)
    min.iter <- as.integer(min.iter)
    min.pairs <- as.integer(min.pairs)
   
    # Checking subset vector and blanking out the unused names.
    subset.row <- .subset2index(subset.row, x, byrow=TRUE)
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
        for (cl in names(pairs)) { 
            message(sprintf("Number of %s pairs: %d", cl, length(pairs[[cl]][[1]])))
        }
    }

    if (.bpNotSharedOrUp(BPPARAM)) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM))
    }

    # Run the allocation algorithm.
    all.scores <- vector('list', length(pairs))
    names(all.scores) <- names(pairs)
    for (cl in names(pairs)) { 
        pcg.state <- .setup_pcg_state(ncol(x))
        pairings <- pairs[[cl]]
        cur.scores <- colBlockApply(x, FUN=.cyclone_scores, niters=iter, miniters=min.iter, 
            minpairs=min.pairs, marker1=pairings$first, marker2=pairings$second, indices=pairings$index,
            seeds=pcg.state$seeds[[1]], streams=pcg.state$streams[[1]], BPPARAM=BPPARAM)
        all.scores[[cl]] <- unlist(cur.scores)
    }

    # Assembling the output.
    scores <- do.call(data.frame, all.scores)
    scores.normalised <- scores/rowSums(scores)

    # Getting the phases.
    phases <- ifelse(scores$G1 >= scores$G2M, "G1", "G2M")
    phases[scores$G1 < 0.5 & scores$G2M < 0.5] <- "S"

    list(phases=phases, scores=scores, normalized.scores=scores.normalised)
}

#' @importFrom DelayedArray currentViewport makeNindexFromArrayViewport
.cyclone_scores <- function(block, ..., seeds, streams) {
    vp <- currentViewport()
    cols <- makeNindexFromArrayViewport(vp, expand.RangeNSBS=TRUE)[[2]]
    if (!is.null(cols)) {
        seeds <- seeds[cols]
        streams <- streams[cols]
    }
    cyclone_scores(block, ..., seeds=seeds, streams=streams)
}

#' @export
#' @rdname cyclone
setMethod("cyclone", "ANY", .cyclone)

#' @export
#' @rdname cyclone
#' @importFrom SummarizedExperiment assay
setMethod("cyclone", "SummarizedExperiment", function(x, ..., assay.type="counts") {
    cyclone(assay(x, i=assay.type), ...)
})
