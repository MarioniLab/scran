#' Combine variance decompositions
#'
#' Combine the results of multiple variance decompositions, usually generated for the same genes across separate batches of cells.
#'
#' @param ... Two or more \linkS4class{DataFrames} of variance modelling results.
#' For \code{combineVar}, these should be produced by \code{\link{modelGeneVar}} or \code{\link{modelGeneVarWithSpikes}}.
#' For \code{combineCV2}, these should be produced by \code{\link{modelGeneCV2}} or \code{\link{modelGeneCV2WithSpikes}}.
#' @param method String specifying how p-values are to be combined, see \code{\link{combinePValues}} for options.
#' @param equiweight Logical scalar indicating whether each result is to be given equal weight in the combined statistics.
#' @param ncells Numeric vector containing the number of cells used to generate each element of \code{...}.
#' Only used if \code{equiweight=FALSE}.
#'
#' @details
#' These functions are designed to merge results from separate calls to \code{\link{modelGeneVar}}, \code{\link{modelGeneCV2}} or related functions, where each result is usually computed for a different batch of cells.
#' Separate variance decompositions are necessary in cases where different concentrations of spike-in have been added to the cells in each batch.
#' This affects the technical mean-variance relationship and precludes the use of a common trend fit.
#' 
#' For \code{combineVar}, the combined mean is computed by averaging the means across all input DataFrames.
#' Similarly, each variance component is combined by taking the average of the corresponding value across inputs.
#' For \code{combineCV2}, the same process is applied but using the geometric mean of the various statistics.
#' This difference reflects the strategy in which the relevant measure of overdispersion is computed - 
#' \code{bio} is computed by subtraction while \code{ratio} is computed by division.
#'
#' If \code{equiweight=FALSE}, each per-batch statistic is weighted by the number of cells used to compute it.
#' The number of cells can be explicitly set using \code{ncells}, and is otherwise assumed to be equal for all batches.
#' No weighting is performed by default, which ensures that all batches contribute equally to the combined statistics.
#' This avoids cases where batches with many cells dominate the output.
#'
#' The \code{\link{combinePValues}} function is used to combine p-values across batches.
#' Only \code{method="z"} will perform any weighting of batches, and only if \code{weights} is set.
#' 
#' @return
#' A DataFrame with the same numeric fields as that produced by \code{\link{modelGeneVar}} or \code{\link{modelGeneCV2}}.
#' Each row corresponds to an input gene.
#' Each field contains the (weighted) arithmetic/geometric mean across all batches except for \code{p.value}, which contains the combined p-value based on \code{method};
#' and \code{FDR}, which contains the adjusted p-value using the BH method.
#' 
#' @seealso
#' \code{\link{modelGeneVar}} and \code{\link{modelGeneCV2}}, for two possible inputs into this function.
#'
#' \code{\link{combinePValues}}, for details on how the p-values are combined.
#' 
#' @author
#' Aaron Lun
#' 
#' @examples
#' data(example.sce)
#' 
#' library(scater)
#' y1 <- example.sce[,1:100] 
#' y1 <- logNormCounts(y1) # normalize separately after subsetting.
#' results1 <- modelGeneVar(y1)
#' 
#' y2 <- example.sce[,1:100 + 100] 
#' y2 <- logNormCounts(y2) # normalize separately after subsetting.
#' results2 <- modelGeneVar(y2)
#' 
#' head(combineVar(results1, results2))
#' head(combineVar(results1, results2, method="simes"))
#' head(combineVar(results1, results2, method="berger"))
#' 
#' @export
combineVar <- function(..., method="fisher", equiweight=TRUE, ncells=NULL) {
    collected <- list(...)
    if (is.null(ncells)) {
        ncells <- rep(10L, length(collected))
    }
    .combine_blocked_statistics(collected, method=method, equiweight=equiweight, ncells=ncells)
}

#' @export
#' @rdname combineVar
combineCV2 <- function(..., method="fisher", equiweight=TRUE, ncells=NULL) {
    collected <- list(...)
    if (is.null(ncells)) {
        ncells <- rep(10L, length(collected))
    }
    .combine_blocked_statistics(collected, method=method, equiweight=equiweight, ncells=ncells,
        geometric=TRUE, fields=c("mean", "total", "trend", "ratio"))
}
