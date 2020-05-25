#' Hypothesis tests with linear models
#'
#' Perform basic hypothesis tests with linear models in an efficient manner.
#'
#' @param x A numeric matrix-like object containing log-expression values for cells (columns) and genes (rows).
#' Alternatively, a \linkS4class{SummarizedExperiment} containing such a matrix.
#' @param design A numeric design matrix with number of rows equal to \code{ncol(x)}.
#' @param contrast A numeric vector or matrix specifying the contrast of interest.
#' This should have length (if vector) or number of rows (if matrix) equal to \code{ncol(x)}.
#' @param ... Further arguments to pass to \code{\link{fitLinearModel}}.
#'
#' @return A \linkS4class{DataFrame} containing test results with one row per row of \code{x}.
#' It contains the estimated values of the contrasted coefficients
#' as well as the p-value and FDR for each gene.
#'
#' @details
#' The null hypothesis is defined by the contrast matrix or vector in the same manner 
#' that is used in the \pkg{limma} and \pkg{edgeR} packages.
#' Briefly, the contrast vector specifies a linear combination of coefficients that sums to zero under the null.
#' For contrast matrices, the joint null consists of the intersection of the nulls defined by each column vector.
#'
#' This function can be considered a more efficient version of \code{\link{lmFit}}
#' that works on a variety of matrix representations (see \code{\link{fitLinearModel}}).
#' It also omits the empirical Bayes shrinkage step,
#' which is acceptable given the large number of residual d.f. in typical single-cell studies.
#'
#' @author Aaron Lun
#'
#' @examples
#' y <- matrix(rnorm(10000), ncol=100)
#'
#' # Example with categorical factors:
#' A <- gl(2, 50)
#' design <- model.matrix(~A)
#' testLinearModel(y, design, contrast=c(0, 1))
#'
#' # Example with continuous variables:
#' u <- runif(100)
#' design <- model.matrix(~u)
#' testLinearModel(y, design, contrast=c(0, 1))
#'
#' # Example with multiple varibales:
#' B <- gl(4, 25)
#' design <- model.matrix(~B)
#' testLinearModel(y, design, contrast=cbind(c(0,1,0,0), c(0,0,1,-1)))
#'
#' @name testLinearModel
NULL

#' @importFrom scuttle fitLinearModel
#' @importFrom limma lmFit classifyTestsF contrasts.fit
#' @importFrom S4Vectors DataFrame
#' @importFrom stats p.adjust pt pf
.test_linear_model <- function(x, design, contrast, ...) {
    full <- fitLinearModel(x, design, ...)

    # Hacking limma to compute our desired statistics.
    lfit <- lmFit(rbind(seq_len(nrow(design))), design)
    lfit$coefficients <- full$coefficients
    lfit$sigma2 <- full$variance
    lfit <- contrasts.fit(lfit, contrast)

    tstat <- lfit$coefficients / outer(sqrt(lfit$sigma2), lfit$stdev.unscaled[1,])

    if (ncol(tstat)==1L) {
        tstat <- drop(tstat)
        pvalue <- pt(abs(tstat), df=full$residual.df, lower.tail=FALSE) * 2
        if (is.null(colnames(lfit$coefficients))) {
            colnames(lfit$coefficients) <- "logFC"
        }
    } else {
        lfit$tstat <- tstat
        fstat <- classifyTestsF(lfit, fstat.only=TRUE)
        pvalue <- pf(fstat, ncol(tstat), full$residual.df, lower.tail = FALSE)
        attributes(pvalue) <- NULL
    }
   
    DataFrame(row.names=rownames(full$coefficients), # account for subsetting.
        lfit$coefficients, 
        p.value=pvalue,
        FDR=p.adjust(pvalue, method="BH"))
}

#' @export
setGeneric("testLinearModel", function(x, ...) standardGeneric("testLinearModel"))

#' @export
setMethod("testLinearModel", "ANY", .test_linear_model)

#' @export
#' @importFrom SummarizedExperiment assay
setMethod("testLinearModel", "SummarizedExperiment", function(x, ..., assay.type="logcounts") {
    .test_linear_model(assay(x, assay.type), ...)
})
