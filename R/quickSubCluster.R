#' Quick and dirty subclustering
#'
#' Performs a quick subclustering for all cells within each group.
#'
#' @param x A matrix of counts or log-normalized expression values (if \code{normalize=FALSE}),
#' where each row corresponds to a gene and each column corresponds to a cell.
#' 
#' Alternatively, a \linkS4class{SummarizedExperiment} or \linkS4class{SingleCellExperiment} object containing such a matrix.
#' @param groups A vector of group assignments for all cells, usually corresponding to cluster identities.
#' @param normalize Logical scalar indicating whether each subset of \code{x} should be log-transformed prior to further analysis.
#' @param prepFUN A function that accepts a single \linkS4class{SingleCellExperiment} object and returns another \linkS4class{SingleCellExperiment} containing any additional elements required for clustering (e.g., PCA results).
#' @param min.ncells An integer scalar specifying the minimum number of cells in a group to be considered for subclustering.
#' @param clusterFUN A function that accepts a single \linkS4class{SingleCellExperiment} object and returns a vector of cluster assignments for each cell in that object.
#' @param format A string to be passed to \code{\link{sprintf}}, specifying how the subclusters should be named with respect to the parent level in \code{groups} and the level returned by \code{clusterFUN}.
#' @param assay.type String or integer scalar specifying the relevant assay.
#' @param ... For the generic, further arguments to pass to specific methods.
#'
#' For the ANY and SummarizedExperiment methods, further arguments to pass to the SingleCellExperiment method.
#'
#' @return
#' A named \linkS4class{List} of \linkS4class{SingleCellExperiment} objects.
#' Each object corresponds to a level of \code{groups} and contains a \code{"subcluster"} column metadata field with the subcluster identities for each cell.
#'
#' The \code{\link{metadata}} of the List also contains \code{index}, a list of integer vectors specifying the cells in \code{x} in each returned SingleCellExperiment object; and \code{subcluster}, a character vector of subcluster identities for all cells in \code{x}.
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
#' % We use denoisePCA+modelGeneVar by default here, because we hope that each parent cluster is reasonably homogeneous.
#' % This allows us to assume that the trend is actually a good estimate of the technical noise.
#' % We don't use the other modelGeneVar*'s to avoid making assumptions about the available of spike-ins, UMI data, etc.
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
#' @importFrom S4Vectors metadata<-
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

    all.sce <- collated <- list()
    by.group <- split(seq_along(groups), groups)

    for (i in names(by.group)) {
        y <- x[,by.group[[i]]]
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
        collated[[i]] <- clusters
    }

    all.clusters <- rep(NA_character_, ncol(x))
    all.clusters[unlist(by.group)] <- unlist(collated)

    output <- as(all.sce, "List")
    metadata(output) <- list(index=by.group, subcluster=all.clusters)
    output
}

############################
# S4 method definitions
############################

# NOTE: normally, I would have dispatch flow "downwards" towards base
# classes or simpler objects. However, in this case, I have implemented
# it to flow upwards because prepFUN and clustFUN expect SCE inputs.
# This means that everything ends up being promoted to an SCE anyway,
# and it's too hard to require users write their functions endomorphically.

#' @export
#' @rdname quickSubCluster
setGeneric("quickSubCluster", function(x, ...) standardGeneric("quickSubCluster"))

#' @export
#' @rdname quickSubCluster
#' @importFrom SingleCellExperiment SingleCellExperiment
setMethod("quickSubCluster", "ANY", function(x, normalize=TRUE, ...)
{
    assays <- list(x)
    assay.type <- if (normalize) "counts" else "logcounts"
    names(assays) <- assay.type 
    .quick_sub_cluster(SingleCellExperiment(assays), normalize=normalize, assay.type=assay.type, ...)
})

#' @export
#' @rdname quickSubCluster
#' @importFrom methods as
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
setMethod("quickSubCluster", "SummarizedExperiment", function(x, ...) {
    .quick_sub_cluster(as(x, "SingleCellExperiment"), ...)
})

#' @export
#' @rdname quickSubCluster
setMethod("quickSubCluster", "SingleCellExperiment", .quick_sub_cluster)
