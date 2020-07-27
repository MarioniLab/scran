#' Assess cluster stability by bootstrapping
#'
#' Generate bootstrap replicates and recluster on them to determine the stability of clusters with respect to sampling noise.
#' This has been migrated to the \pkg{bluster} package as \code{\link{bootstrapStability}}.
#'
#' @param ... Arguments to pass to \code{\link{bootstrapStability}}.
#'
#' @return 
#' See \code{\link{bootstrapStability}} for details.
#' 
#' @author
#' Aaron Lun
#' 
#' @seealso
#' \code{\link{bootstrapStability}}, the new name for this function.
#' 
#' @export
#' @importFrom bluster bootstrapStability
bootstrapCluster <- function(...) {
    .Deprecated(new="bluster::bootstrapStability")
    bootstrapStability(...)
}
