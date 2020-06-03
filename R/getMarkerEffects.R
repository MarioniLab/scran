#' Get marker effect sizes
#'
#' Utility function to extract the marker effect sizes as a matrix from the output of \code{\link{findMarkers}}.
#'
#' @param x A \linkS4class{DataFrame} containing marker statistics for a given group/cluster,
#' usually one element of the List returned by \code{\link{findMarkers}}.
#' @param prefix String containing the prefix for the columns containing the effect size.
#' @param strip Logical scalar indicating whether the prefix should be removed from the output column names.
#' @param remove.na.col Logical scalar indicating whether to remove columns containing any \code{NA}s.
#'
#' @details
#' Setting \code{remove.na.col=TRUE} may be desirable in applications involving blocked comparisons,
#' where some pairwise comparisons are not possible if the relevant levels occur in different blocks.
#' In such cases, the resulting column is filled with \code{NA}s that may interfere with downstream steps like clustering.
#' 
#' @return A numeric matrix containing the effect sizes for the comparison to every other group/cluster.
#' 
#' @author Aaron Lun
#'
#' @examples
#' library(scuttle)
#' sce <- mockSCE()
#' sce <- logNormCounts(sce)
#'                                                                      
#' kout <- kmeans(t(logcounts(sce)), centers=4) 
#' out <- findMarkers(sce, groups=kout$cluster)
#'
#' eff1 <- getMarkerEffects(out[[1]])
#' str(eff1)
#'
#' @seealso
#' \code{\link{findMarkers}} and \code{\link{combineMarkers}}, to generate the DataFrames.
#'
#' @export
#' @importFrom DelayedMatrixStats colAnyNAs
getMarkerEffects <- function(x, prefix="logFC", strip=TRUE, remove.na.col=FALSE) {
    regex <- paste0("^", prefix, "\\.")
    i <- grep(regex, colnames(x))
    out <- as.matrix(x[,i])

    if (strip) {
        colnames(out) <- sub(regex, "", colnames(out))
    }
    if (remove.na.col) {
        out <- out[,!colAnyNAs(out),drop=FALSE]
    }

    out
}
