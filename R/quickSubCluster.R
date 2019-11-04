#' Quick and dirty subclustering
#'
#' Performs a quick subclustering for all cells within each group.
#'
#' @param x A matrix of counts or log-normalized expression values (if \code{normalize=FALSE}).
#' Alternatively, a \linkS4class{SingleCellExperiment} object containing such a matrix.
#' @param groups A vector of group assignments for all cells, usually corresponding to cluster identities.
#' @param normalize Logical scalar indicating whether each subset of \code{x} should be log-transformed prior to further analysis.
#' @param prepFUN A function that accepts a single \linkS4class{SingleCellExperiment} object and returns another \linkS4class{SingleCellExperiment} containing any additional elements required for clustering (e.g., PCA results).
#' @param min.ncells An integer scalar specifying the minimum number of cells in a group to be considered for subclustering.
#' @param clusterFUN A function that accepts a single \linkS4class{SingleCellExperiment} object and returns a vector of cluster assignments for each cell in that object.
#' @param format A string to be passed to \code{\link{sprintf}}, specifying how the subclusters should be named with respect to the parent level in \code{groups} and the level returned by \code{clusterFUN}.
#' @param assay.type String or integer scalar specifying the relevant assay.
#' @param ... Further arguments to pass to specific methods.
#'
#' @return
#' A named \linkS4class{List} of \linkS4class{SingleCellExperiment} objects.
#' Each object corresponds to a level of \code{groups} and contains a \code{"subcluster"} column metadata field with the subcluster identities for each cell.
#'
#' @details
#' \code{quickSubCluster} is a simple convenience function that loops over all levels of \code{groups} to perform subclustering.
#' It subsets \code{x} to retain all cells in one level and then runs \code{prepFUN} and \code{clusterFUN} to cluster them.
#' Levels with fewer than \code{min.ncells} are not subclustered and have \code{"subcluster"} set to the name of the level.
#'
#' The distinction between \code{prepFUN} and \code{clusterFUN} is that the former's calculations are preserved in the output.
#' For example, we would put the PCA in \code{prepFUN} so that the PCs are returned in the \code{\link{reducedDims}} for later use.
#' In contrast, \code{clusterFUN} is only used to obtain the subcluster assignments so any intermediate objects are lost.
#' 
#' By default, \code{prepFUN} will run \code{\link{modelGeneVar}}, take the top 10% of genes with large biological components with \code{\link{getTopHVGs}}, and then run \code{\link{denoisePCA}} to perform the PCA.
#' \code{clusterFUN} will then perform graph-based clustering with \code{\link{buildSNNGraph}} and \code{\link{cluster_walktrap}}.
#' Either or both of these functions can be replaced with custom functions.
#'
#' The default behavior of this function is the same as running \code{\link{quickCluster}} on each subset with default parameters except for \code{min.size=0}.
#'
#' @author Aaron Lun
#' 
#' @examples
#' library(scater)
#' sce <- mockSCE(ncells=200)
#' 
#' # Lowering min.size for this small demo:
#' clusters <- quickCluster(sce, min.size=50)
#'
#' # Getting subclusters:
#' out <- quickSubCluster(sce, clusters)
#'
#' # Defining custom prep functions:
#' out2 <- quickSubCluster(sce, clusters, 
#'     prepFUN=function(x) {
#'         dec <- modelGeneVarWithSpikes(x, "Spikes")
#'         top <- getTopHVGs(dec, prop=0.2)
#'         runPCA(x, subset_row=top, ncomponents=25)
#'     }
#' )
#' 
#' # Defining custom cluster functions:
#' out3 <- quickSubCluster(sce, clusters, 
#'     clusterFUN=function(x) {
#'         kmeans(reducedDim(x, "PCA"), sqrt(ncol(x)))$cluster
#'     }
#' )
#'
#' @seealso
#' \code{\link{quickCluster}}, for a related function to quickly obtain clusters.
#'
#' @name quickSubCluster
NULL

#' @importFrom scater logNormCounts
#' @importFrom igraph cluster_walktrap
#' @importFrom BiocSingular bsparam
#' @importClassesFrom S4Vectors List
.quick_sub_cluster <- function(x, groups, normalize=TRUE, 
    prepFUN=NULL, min.ncells=50, clusterFUN=NULL, 
    format="%s.%s", assay.type="counts") 
{
    if (normalize) {
        alt.assay <- "logcounts"
    } else {
        alt.assay <- assay.type
    }

    if (is.null(prepFUN)) {
        prepFUN <- function(x) {
            # Putting this check here to avoid skipping
            # user-specified modifications.
            if (ncol(x) < 2L) {
                return(x)
            }

            # For consistency with quickCluster().
            dec <- modelGeneVar(x, assay.type=alt.assay)
            top <- getTopHVGs(dec, n=500, prop=0.1)
            denoisePCA(x, dec, subset.row=top, assay.type=alt.assay)
        }
    }
    if (is.null(clusterFUN)) {
        clusterFUN <- function(x) {
            g.trans <- buildSNNGraph(x, use.dimred="PCA")
            cluster_walktrap(g.trans)$membership
        }
    }

    all.sce <- list()
    for (i in as.character(sort(unique(groups)))) {
        y <- x[,groups==i]
        if (normalize) {
            y <- logNormCounts(y, exprs_values=assay.type)
        }

        y <- prepFUN(y)
        if (ncol(y) >= min.ncells) {
            clusters <- clusterFUN(y)
            clusters <- sprintf(format, i, clusters)
        } else {
            clusters <- rep(i, ncol(y)) 
        }

        y$subcluster <- clusters
        all.sce[[i]] <- y
    }

    as(all.sce, "List")
}

############################
# S4 method definitions
############################

#' @export
#' @rdname quickSubCluster
setGeneric("quickSubCluster", function(x, ...) standardGeneric("quickSubCluster"))

#' @export
#' @rdname quickSubCluster
#' @importFrom SingleCellExperiment SingleCellExperiment
setMethod("quickSubCluster", "ANY", function(x, groups, normalize=TRUE, 
    prepFUN=NULL, min.ncells=50, clusterFUN=NULL, format="%s.%s") 
{
    assays <- list(x)
    assay.type <- if (normalize) "counts" else "logcounts"
    names(assays) <- assay.type
    .quick_sub_cluster(SingleCellExperiment(assays), groups, 
        normalize=normalize, prepFUN=prepFUN, min.ncells=min.ncells, 
        clusterFUN=clusterFUN, format=format, assay.type=assay.type)
})

#' @export
#' @rdname quickSubCluster
#' @importFrom SummarizedExperiment assay
setMethod("quickSubCluster", "SingleCellExperiment", .quick_sub_cluster)