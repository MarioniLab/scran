#' Normalization by deconvolution
#'
#' Scaling normalization of single-cell RNA-seq data by deconvolving size factors from cell pools.
#' These functions have been moved to the \pkg{scuttle} package and are just retained here for compatibility.
#'
#' @param ... Further arguments to pass to \code{\link{pooledSizeFactors}} or \code{\link{computePooledFactors}}.
#'
#' @return
#' For \code{calculateSumFactors}, a numeric vector of size factors returned by \code{\link{pooledSizeFactors}}.
#'
#' For \code{computeSumFactors}, a SingleCellExperiment containing the size factors in its \code{\link{sizeFactors}},
#' as returned by \code{\link{computePooledFactors}}.
#'
#' @author Aaron Lun
#' @export
#' @importFrom scuttle computePooledFactors
computeSumFactors <- function(...) {
    computePooledFactors(...)
}

#' @export
#' @importFrom scuttle pooledSizeFactors
#' @rdname computeSumFactors 
calculateSumFactors <- function(...) {
    pooledSizeFactors(...)
}
