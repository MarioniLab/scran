#' Compute the cluster-wise modularity
#' 
#' Calculate the modularity of each cluster from a graph, based on a null model of random connections between nodes.
#' This function has been moved to the \pkg{bluster} package as \code{\link{pairwiseModularity}}.
#'
#' @param ... Arguments to pass to \code{\link{pairwiseModularity}}.
#' 
#' @return
#' See \code{?\link{pairwiseModularity}} for more details.
#'
#' @author
#' Aaron Lun
#' 
#' @seealso
#' \code{\link{pairwiseModularity}}, the new name for this function.
#' 
#' @export
#' @importFrom bluster pairwiseModularity
clusterModularity <- function(...) {
    .Deprecated(new="bluster::pairwiseModularity")
    pairwiseModularity(...)
}
