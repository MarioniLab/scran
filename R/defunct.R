#' Defunct functions
#'
#' Functions that have passed on to the function afterlife.
#' Their successors are also listed.
#'
#' @param ... Ignored arguments.
#'
#' @section Variance modelling:
#' \code{trendVar}, \code{decomposeVar} and \code{testVar} are succeeded by a suite of funtions related to \code{\link{modelGeneVar}} and \code{\link{fitTrendVar}}.
#'
#' \code{improvedCV2} and \code{technicalCV2} are succeeded by \code{\link{modelGeneCV2}} and \code{\link{fitTrendCV2}}.
#'
#' \code{makeTechTrend} is succeeded by \code{\link{modelGeneVarByPoisson}}.
#'
#' \code{multiBlockVar} is succeeded by the \code{block} argument in many of the modelling functions, and \code{multiBlockNorm} is no longer necessary.
#'
#' @section Clustering-related functions:
#' \code{bootstrapCluster} has been moved over to the \pkg{bluster} package, as the \code{\link{bootstrapStability}} function.
#'
#' \code{neighborsToSNNGraph} and \code{neighborsToKNNGraph} have been moved over to the \pkg{bluster} package.
#'
#' \code{clusterModularity} has been moved over to the \pkg{bluster} package, as the \code{\link{pairwiseModularity}} function.
#'
#' \code{clusterPurity} has been moved over to the \pkg{bluster} package, as the \code{\link{neighborPurity}} function.
#'
#' \code{clusterSNNGraph} and \code{clusterKNNGraph} have been replaced by \code{\link{clusterRows}} with \linkS4class{NNGraphParam} or \linkS4class{TwoStepParam} from the \pkg{bluster} package.
#'
#' \code{coassignProb} and \code{clusterRand} have been replaced by \code{\link{pairwiseRand}} from the \pkg{bluster} package.
#'
#' @section Pseudotime-related functions:
#' \code{createClusterMST}, \code{quickPseudotime} and \code{testPseudotime} have been moved over to the \pkg{TSCAN} package.
#'
#' \code{connectClusterMST} has been moved over to the \pkg{TSCAN} package, as the \code{reportEdges} function.
#'
#' \code{orderClusterMST} has been moved over to the \pkg{TSCAN} package, as the \code{orderCells} function.
#'
#' @section Doublet-related functions:
#' \code{doubletCells} has been moved over to the \pkg{scDblFinder} package, as the \code{computeDoubletDensity} function.
#'
#' \code{doubletCluster} has been moved over to the \pkg{scDblFinder} package, as the \code{findDoubletClusters} function.
#'
#' \code{doubletRecovery} has been moved over to the \pkg{scDblFinder} package, as the \code{recoverDoublets} function.
#'
#' @section Other functions:
#' \code{overlapExprs} is succeeded by \code{\link{findMarkers}} with \code{test.type="wilcox"}.
#'
#' \code{parallelPCA} has been moved over to the \pkg{PCAtools} package.
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

#' @export
#' @rdname defunct
bootstrapCluster <- function(...) {
    .Defunct("bluster::bootstrapStability")
}

#' @export
#' @rdname defunct
clusterModularity <- function(...) {
    .Defunct("bluster::pairwiseModularity")
}

#' @export
#' @rdname defunct
clusterPurity <- function(...) {
    .Defunct("bluster::neighborPurity")
}

#' @export
#' @rdname defunct
clusterKNNGraph <- function(...) {
    .Defunct(msg="'clusterKNNGraph' is deprecated.\nUse 'bluster::clusterRows' with 'NNGraphParam' or 'TwoStepParam' instead.")
}

#' @export
#' @rdname defunct
clusterSNNGraph <- function(...) {
    .Defunct(msg="'clusterSNNGraph' is deprecated.\nUse 'bluster::clusterRows' with 'NNGraphParam' or 'TwoStepParam' instead.")
}

#' @export
#' @rdname defunct
coassignProb <- function(...) {
    .Defunct("bluster::pairwiseRand")
}

#' @export
#' @rdname defunct
createClusterMST <- function(...) {
    .Defunct("TSCAN::createClusterMST")
}

#' @export
#' @rdname defunct
connectClusterMST <- function(...) {
    .Defunct("TSCAN::reportEdges")
}

#' @export
#' @rdname defunct
orderClusterMST <- function(...) {
    .Defunct("TSCAN::mapCellsToEdges")
}

#' @export
#' @rdname defunct
quickPseudotime <- function(...) {
    .Defunct("TSCAN::quickPseudotime")
}

#' @export
#' @rdname defunct
testPseudotime <- function(...) {
    .Defunct("TSCAN::testPseudotime")
}

#' @export
#' @rdname defunct
doubletCells <- function(...) {
    .Defunct("scDblFinder::computeDoubletDensity")
}

#' @export
#' @rdname defunct
doubletCluster <- function(...) {
    .Defunct("scDblFinder::findDoubletClusters")
}

#' @export
#' @rdname defunct
doubletRecovery <- function(...) {
    .Defunct("scDblFinder::recoverDoublets")
}

