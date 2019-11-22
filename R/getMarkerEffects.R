#' Get marker effect sizes
#'
#' Utility function to extract the marker effect sizes as a matrix from the output of \code{\link{findMarkers}}.
#'
#' @param x A \linkS4class{DataFrame} containing marker statistics for a given group/cluster,
#' usually one element of the List returned by \code{\link{findMarkers}}.
#' @param prefix String containing the prefix for the columns containing the effect size.
#' @param strip Logical scalar indicating whether the prefix should be removed from the output column names.
#'
#' @return A numeric matrix containing the effect sizes for the comparison to every other group/cluster.
#' 
#' @author Aaron Lun
#'
#' @examples
#' library(scater)
#' sce <- mockSCE()
#' sce <- logNormCounts(sce)
#'                                                                      
#' kout <- kmeans(t(logcounts(sce)), centers=4) 
#' out <- findMarkers(sce, groups=kout$cluster)
#'
#' eff1 <- getMarkerEffects(out[[1]])
#' str(eff1)
#' @seealso
#' \code{\link{findMarkers}} and \code{\link{combineMarkers}}, to generate the DataFrames.
#'
#' @export
getMarkerEffects <- function(x, prefix="logFC", strip=TRUE) {
    regex <- paste0("^", prefix, "\\.")
    i <- grep(regex, colnames(x))
    out <- as.matrix(x[,i])

    if (strip) {
        colnames(out) <- sub(regex, "", colnames(out))
    }
    out
}
