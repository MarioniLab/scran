#' Combine variance decompositions
#'
#' Combine the results of multiple variance decompositions, usually generated for the same genes across separate batches of cells.
#'
#' @param ... Two or more \linkS4class{DataFrame}s of variance modelling results.
#' For \code{combineVar}, these should be produced by \code{\link{modelGeneVar}} or \code{\link{modelGeneVarWithSpikes}}.
#' For \code{combineCV2}, these should be produced by \code{\link{modelGeneCV2}} or \code{\link{modelGeneCV2WithSpikes}}.
#'
#' Alternatively, one or more lists of DataFrames containing variance modelling results.
#' Mixed inputs are also acceptable, e.g., lists of DataFrames alongside the DataFrames themselves.
#' @param method String specifying how p-values are to be combined, see \code{\link{combineParallelPValues}} for options.
#' @param pval.field A string specifying the column name of each element of \code{...} that contains the p-value.
#' @param other.fields A character vector specifying the fields containing other statistics to combine.
#' @param equiweight Logical scalar indicating whether each result is to be given equal weight in the combined statistics.
#' @param ncells Numeric vector containing the number of cells used to generate each element of \code{...}.
#' Only used if \code{equiweight=FALSE}.
#'
#' @details
#' These functions are designed to merge results from separate calls to \code{\link{modelGeneVar}}, \code{\link{modelGeneCV2}} or related functions, where each result is usually computed for a different batch of cells.
#' Separate variance decompositions are necessary in cases where the mean-variance relationships vary across batches (e.g., different concentrations of spike-in have been added to the cells in each batch), which precludes the use of a common trend fit.
#' By combining these results into a single set of statistics, we can apply standard strategies for feature selection in multi-batch integrated analyses. 
#' 
#' By default, statistics in \code{other.fields} contain all common non-numeric fields that are not \code{pval.field} or \code{"FDR"}.
#' This usually includes \code{"mean"}, \code{"total"}, \code{"bio"} (for \code{combineVar}) or \code{"ratio"} (for \code{combineCV2}).
#' \itemize{
#' \item For \code{combineVar}, statistics are combined by averaging them across all input DataFrames.
#' \item For \code{combineCV2}, statistics are combined by taking the geometric mean across all inputs.
#' }
#' This difference between functions reflects the method by which the relevant measure of overdispersion is computed.
#' For example, \code{"bio"} is computed by subtraction, so taking the average \code{bio} remains consistent with subtraction of the total and technical averages.
#' Similarly, \code{"ratio"} is computed by division, so the combined \code{ratio} is consistent with division of the geometric means of the total and trend values.
#'
#' If \code{equiweight=FALSE}, each per-batch statistic is weighted by the number of cells used to compute it.
#' The number of cells can be explicitly set using \code{ncells}, and is otherwise assumed to be equal for all batches.
#' No weighting is performed by default, which ensures that all batches contribute equally to the combined statistics and avoids situations where batches with many cells dominate the output.
#'
#' The \code{\link{combineParallelPValues}} function is used to combine p-values across batches.
#' The default is to use Fisher's method, which will achieve a low p-value if a gene is highly variable in any batch.
#' Only \code{method="stouffer"} will perform any weighting of batches, and only if \code{weights} is set.
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
#' \code{\link{combineParallelPValues}}, for details on how the p-values are combined.
#' 
#' @author
#' Aaron Lun
#' 
#' @examples
#' library(scuttle)
#' sce <- mockSCE()
#' 
#' y1 <- sce[,1:100] 
#' y1 <- logNormCounts(y1) # normalize separately after subsetting.
#' results1 <- modelGeneVar(y1)
#' 
#' y2 <- sce[,1:100 + 100] 
#' y2 <- logNormCounts(y2) # normalize separately after subsetting.
#' results2 <- modelGeneVar(y2)
#' 
#' head(combineVar(results1, results2))
#' head(combineVar(results1, results2, method="simes"))
#' head(combineVar(results1, results2, method="berger"))
#' 
#' @export
#' @importFrom scuttle .unpackLists
combineVar <- function(..., method="fisher", pval.field="p.value", other.fields=NULL, equiweight=TRUE, ncells=NULL) {
    collected <- .unpackLists(...)
    if (is.null(ncells)) {
        # Any arbitrary value will do here, 
        # as long as it's >=2 so that the block isn't ignored.
        ncells <- rep(10L, length(collected))
    }
    if (is.null(other.fields)) {
        other.fields <- .find_other_fields(collected, c(pval.field, "FDR"))
    }
    .combine_blocked_statistics(collected, method=method, equiweight=equiweight, ncells=ncells,
        geometric=FALSE, fields=other.fields, pval=pval.field)                                
}

#' @export
#' @rdname combineVar
#' @importFrom scuttle .unpackLists
combineCV2 <- function(..., method="fisher", pval.field="p.value", other.fields=NULL, equiweight=TRUE, ncells=NULL) {
    collected <- .unpackLists(...)
    if (is.null(ncells)) {
        ncells <- rep(10L, length(collected))
    }
    if (is.null(other.fields)) {
        other.fields <- .find_other_fields(collected, c(pval.field, "FDR"))
    }
    .combine_blocked_statistics(collected, method=method, equiweight=equiweight, ncells=ncells,
        geometric=TRUE, fields=other.fields, pval=pval.field)
}

.find_other_fields <- function(collected, exclude) {
    all.numerics <- list()
    for (i in seq_along(collected)) {
        x <- collected[[i]]
        all.numerics[[i]] <- colnames(x)[vapply(x, is.numeric, TRUE)]
    }
    setdiff(Reduce(intersect, all.numerics), exclude)
}
