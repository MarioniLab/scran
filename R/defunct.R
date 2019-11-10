#' Defunct functions
#'
#' Functions that have passed on to the function afterlife.
#' Their successors are also listed.
#'
#' @param ... Ignored arguments.
#'
#' @details
#' \code{trendVar}, \code{decomposeVar} and \code{testVar} are succeeded by a suite of funtions related to \code{\link{modelGeneVar}} and \code{\link{fitTrendVar}}.
#'
#' \code{improvedCV2} and \code{technicalCV2} are succeeded by \code{\link{modelGeneCV2}} and \code{\link{fitTrendCV2}}.
#'
#' \code{makeTechTrend} is succeeded by \code{\link{modelGeneVarByPoisson}}.
#'
#' \code{multiBlockVar} is succeeded by the \code{block} argument in many of the modelling functions, and \code{multiBlockNorm} is no longer necessary.
#'
#' \code{\link{overlapExprs}} is succeeded by \code{\link{findMarkers}} with \code{test.type="wilcox"}.
#'
#' \code{\link{parallelPCA}} has been moved over to the \pkg{PCAtools} package.
#'
#' @return All functions error out with a defunct message pointing towards its descendent (if available).
#'
#' @author Aaron Lun
#'
#' @examples
#' try(trendVar())
#' @name defunct
NULL

#' @export
#' @rdname defunct
trendVar <- function(...) {
    .Defunct("fitTrendVar")
}

#' @export
#' @rdname defunct
decomposeVar <- function(...) {
    .Defunct("modelGeneVar")
}

#' @export
#' @rdname defunct
testVar <- function(...) {
    .Defunct("modelGeneVar")
}

#' @export
#' @rdname defunct
improvedCV2 <- function(...) {
    .Defunct("modelGeneCV2")
}

#' @export
#' @rdname defunct
technicalCV2 <- function(...) {
    .Defunct("modelGeneCV2")
}

#' @export
#' @rdname defunct
makeTechTrend <- function(...) {
    .Defunct("modelGeneVarByPoisson")
}

#' @export
#' @rdname defunct
multiBlockVar <- function(...) {
    .Defunct()
}

#' @export
#' @rdname defunct
multiBlockNorm <- function(...) {
    .Defunct()
}

#' @export
#' @rdname defunct
overlapExprs <- function(...) {
    .Defunct("findMarkers")
}

#' @export
#' @rdname defunct
parallelPCA <- function(...) {
    .Defunct(msg="'parallelPCA' is defunct.\nUse 'PCAtools::parallelPCA' instead")
}
