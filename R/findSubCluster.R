#' Find subclusters under one or more cluster
#'
#' Performs a quick subclustering for all cells under one or more cluster.
#'
#' @param x A matrix of counts or log-normalized expression values (if
#' \code{normalize=FALSE}), where each row corresponds to a gene and each column
#' corresponds to a cell.
#'
#' Alternatively, a \linkS4class{SummarizedExperiment} or
#' \linkS4class{SingleCellExperiment} object containing such a matrix.
#' @param groups A vector of group assignments for all cells, usually
#' corresponding to cluster identities. If NULL, will determined by
#' \code{colLabels(x)}
#' @param clusters the cluster to be sub-clustered
#' @param normalize Logical scalar indicating whether each subset of \code{x}
#' should be log-transformed prior to further analysis.
#' @param prepFUN A function that accepts a single
#' \linkS4class{SingleCellExperiment} object and returns another
#' \linkS4class{SingleCellExperiment} containing any additional elements
#' required for clustering (e.g., PCA results).
#' @param clusterFUN A function that accepts a single
#' \linkS4class{SingleCellExperiment} object and returns a vector of cluster
#' assignments for each cell in that object.
#' @param BLUSPARAM A \linkS4class{BlusterParam} object that is used to specify
#' the clustering via \code{\link{clusterRows}}.  Only used when
#' \code{clusterFUN=NULL}.
#' @param format A string to be passed to \code{\link{sprintf}}, specifying how
#' the subclusters should be named with respect to the parent level in
#' \code{groups} and the level returned by \code{clusterFUN}.
#' @param assay.type String or integer scalar specifying the relevant assay.
#' @param ... For the generic, further arguments to pass to specific methods.
#'
#' For the ANY and SummarizedExperiment methods, further arguments to pass to
#' the SingleCellExperiment method.
#'
#' @return a factor of length equal to ncol(x) containing cluster assignments
#' for each column of x.
#'
#' @name findSubCluster
NULL

.find_sub_cluster <- function(x, groups = NULL, clusters = NULL, normalize = TRUE, prepFUN = NULL, clusterFUN = NULL, BLUSPARAM = NNGraphParam(), format = "%s.%s", assay.type = NULL) {
    if (is.null(groups)) {
        groups <- SingleCellExperiment::colLabels(x, onAbsence = "error")
    } else if (!identical(length(groups), ncol(x))) {
        stop("the length of groups should equal to ncol(x)")
    }
    # coerce groups and clusters into character vector
    # in case of error matching and indexing
    # keep the original factor levels
    subcluster_levels <- levels(groups)
    groups <- as.character(groups)
    if (is.null(subcluster_levels)) {
        subcluster_levels <- unique(groups)
    }
    if (is.null(clusters)) {
        clusters <- unique(groups)
    } else if (!all(clusters %in% unique(groups))) {
        stop("all clusters specified should be in the groups")
    } else {
        clusters <- as.character(clusters)
    }
    subclusters <- groups
    for (i in clusters) {
        idx <- i == groups
        sce_obj <- x[, idx]
        subcluster_labels <- subcluster_internal(
            sce_obj,
            normalize = normalize, assay.type = assay.type,
            prepFUN = prepFUN, clusterFUN = clusterFUN,
            BLUSPARAM = BLUSPARAM
        )
        subcluster_labels <- sprintf(
            fmt = format,
            subclusters[idx],
            subcluster_labels
        )
        subclusters[idx] <- subcluster_labels
        # insert the new labels into the original levels
        subcluster_levels_idx <- which(i == subcluster_levels)
        subcluster_levels <- append(
            subcluster_levels, 
            sort(unique(subcluster_labels)),
            after = subcluster_levels_idx
        )
        subcluster_levels <- subcluster_levels[-subcluster_levels_idx]
    }
    factor(subclusters, subcluster_levels)
}

# the internal function to implement subcluster
# return a factor of length equal to ncol(x)
subcluster_internal <- function(x, normalize = TRUE, assay.type = NULL, prepFUN = NULL, clusterFUN = NULL, BLUSPARAM = NNGraphParam()) {
    # nolint
    if (is.null(assay.type)) {
        if (normalize) {
            assay.type <- "counts"
        } else {
            assay.type <- "logcounts"
        }
    }
    if (is.null(prepFUN)) {
        prepFUN <- function(x) {
            if (ncol(x) < 2L) {
                return(x)
            }
            dec <- modelGeneVar(x, assay.type = assay.type)
            top <- getTopHVGs(dec, n = 500, prop = 0.1)
            denoisePCA(x, dec, subset.row = top, assay.type = assay.type)
        }
    }
    if (is.null(clusterFUN)) {
        clusterFUN <- function(x) {
            clusterRows(reducedDim(x, "PCA"), BLUSPARAM)
        }
    }
    if (normalize) {
        x <- logNormCounts(x, exprs_values = assay.type)
    }
    x <- prepFUN(x)
    as.character(clusterFUN(x))
}

#' @export
#' @rdname findSubCluster
setGeneric("findSubCluster", function(x, ...) standardGeneric("findSubCluster"))

#' @export
#' @rdname findSubCluster
#' @importFrom SingleCellExperiment SingleCellExperiment
setMethod("findSubCluster", "ANY", function(x, normalize = TRUE, ...) {
    assays <- list(x)
    assay.type <- if (normalize) "counts" else "logcounts"
    names(assays) <- assay.type
    .find_sub_cluster(
        SingleCellExperiment(assays),
        normalize = normalize,
        assay.type = assay.type,
        ...
    )
})

#' @export
#' @rdname findSubCluster
#' @importFrom methods as
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
setMethod("findSubCluster", "SummarizedExperiment", function(x, ...) {
    .find_sub_cluster(as(x, "SingleCellExperiment"), ...)
})

#' @export
#' @rdname findSubCluster
setMethod("findSubCluster", "SingleCellExperiment", .find_sub_cluster)
