#' Combine blockwise statistics
#'
#' Combine DataFrames of statistics computed separately for each block.
#' This usually refers to feature-level statistics and sample-level blocks.
#'
#' @param blocks A list of \linkS4class{DataFrame}s containing blockwise statistics.
#' These should have the same number of rows and the same set of columns.
#' @param ave.fields Character vector specifying the columns of \code{blocks} to be averaged.
#' The value of each column is averaged across blocks, potentially in a weighted manner.
#' @param pval.field String specifying the column of \code{blocks} containing the p-value.
#' This is combined using \code{\link{combineParallelPValues}}.
#' @param method String specifying how p-values should be combined, see \code{?\link{combineParallelPValues}}.
#' @param geometric Logical scalar indicating whether the geometric mean should be computed when averaging \code{ave.fields}.
#' @param equiweight Logical scalar indicating whether each block should be given equal weight.
#' @param weights Numeric vector of length equal to \code{blocks}, containing the weight for each block.
#' Only used if \code{equiweight=TRUE}.
#' @param valid Logical vector indicating whether each block is valid.
#' Invalid blocks are still stored in the \code{per.block} output but are not used to compute the combined statistics.
#'
#' @return A \linkS4class{DataFrame} containing all fields in \code{ave.fields} and the p-values,
#' where each column is created by combining the corresponding block-specific columns.
#' A \code{per.block} column is also reported, containing a DataFrame of the DataFrames of blockwise statistics.
#'
#' @author Aaron Lun
#'
#' @seealso
#' This function is used in \code{\link{modelGeneVar}} and friends, \code{\link{combineVar}} and \code{\link{testLinearModel}}.
#'
#' @examples
#' library(scuttle)
#' sce <- mockSCE()
#' 
#' y1 <- sce[,1:100] 
#' y1 <- logNormCounts(y1) # normalize separately after subsetting.
#' results1 <- modelGeneVar(y1)
#' 
#' y2 <- sce[,1:100 + 100] 
#' y2 <- logNormCounts(y2) # normalize separately after subsetting.
#' results2 <- modelGeneVar(y2)
#'
#' # A manual implementation of combineVar:
#' combineBlocks(list(results1, results2), 
#'     ave.fields=c("mean", "total", "bio", "tech"),
#'     pval.field='p.value', 
#'     method='fisher',
#'     geometric=FALSE,
#'     equiweight=TRUE,
#'     weights=NULL,
#'     valid=c(TRUE, TRUE))
#' 
#' @export
#' @importFrom stats p.adjust
#' @importFrom S4Vectors DataFrame I
#' @importFrom metapod combineParallelPValues
combineBlocks <- function(blocks, ave.fields, pval.field, method, geometric, equiweight, weights, valid) {
    if (length(blocks)==1L) {
        return(blocks[[1]])
    }

    rn <- unique(lapply(blocks, rownames))
    if (length(rn)!=1L) {
        stop("gene identities should be the same")
    }

    if (equiweight) {
        weights <- rep(1, length(blocks))
    } else if (is.null(weights)) {
        stop("'weights' must be specified if 'equiweight=FALSE'")
    }

    original <- blocks
    if (length(unique(vapply(original, nrow, 0L)))!=1L) {
        stop("not all 'blocks' have the same number of rows")
    }

    if (!any(valid)) {
        stop("no entry of 'blocks' has positive weights")
    }
    blocks <- blocks[valid]
    weights <- weights[valid]

    combined <- list()
    for (i in ave.fields) {
        extracted <- lapply(blocks, "[[", i=i)

        if (geometric) {
            extracted <- lapply(extracted, log)
        }
        extracted <- mapply("*", extracted, weights, SIMPLIFY=FALSE, USE.NAMES=FALSE)
        averaged <- Reduce("+", extracted)/sum(weights)
        if (geometric) {
            averaged <- exp(averaged)            
        }
        combined[[i]] <- averaged 
    }

    extracted <- lapply(blocks, "[[", i=pval.field)

    if (method=="z") {
        .Deprecated(old='method="z"', new='method="stouffer"')
        method <- "stouffer"
    } else if (method=="holm-middle") {
        .Deprecated(old='method="holm-middle"', new='method="holm-min"')
        method <- "holm-min"
    }
    combined$p.value <- combineParallelPValues(extracted, method=method, weights=weights)$p.value
    combined$FDR <- p.adjust(combined$p.value, method="BH")

    output <- DataFrame(combined, row.names=rn[[1]])
    output$per.block <- do.call(DataFrame, c(lapply(original, I), list(check.names=FALSE)))

    output
}
