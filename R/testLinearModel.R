#' Hypothesis tests with linear models
#'
#' Perform basic hypothesis tests with linear models in an efficient manner.
#'
#' @param x A numeric matrix-like object containing log-expression values for cells (columns) and genes (rows).
#' Alternatively, a \linkS4class{SummarizedExperiment} containing such a matrix.
#' @param design A numeric design matrix with number of rows equal to \code{ncol(x)}.
#' @param coefs An integer vector specifying the coefficients to drop to form the null model.
#' Only used if \code{contrasts} is not specified.
#' @param contrasts A numeric vector or matrix specifying the contrast of interest.
#' This should have length (if vector) or number of rows (if matrix) equal to \code{ncol(x)}.
#' @param ... For the generic, further arguments to pass to specific methods.
#'
#' For the SummarizedExperiment method, further arguments to pass to the ANY method.
#' @inheritParams modelGeneVar
#'
#' @return A \linkS4class{DataFrame} containing test results with one row per row of \code{x}.
#' It contains the estimated values of the contrasted coefficients
#' as well as the p-value and FDR for each gene.
#'
#' @details
#' This function can be considered a more efficient version of \code{\link{lmFit}}
#' that works on a variety of matrix representations (see \code{\link{fitLinearModel}}).
#' It also omits the empirical Bayes shrinkage step,
#' which is acceptable given the large number of residual d.f. in typical single-cell studies.
#'
#' If \code{contrasts} is specified, the null hypothesis is defined by the contrast matrix or vector in the same manner 
#' that is used in the \pkg{limma} and \pkg{edgeR} packages.
#' Briefly, the contrast vector specifies a linear combination of coefficients that sums to zero under the null.
#' For contrast matrices, the joint null consists of the intersection of the nulls defined by each column vector.
#'
#' Otherwise, if only \code{coefs} is specified, 
#' the null model is formed by simply dropping all of the specified coefficients from \code{design}.
#'
#' If \code{block} is specified, a linear model is fitted separately to the cells in each level.
#' The results are combined across levels by averaging coefficients and combining p-values with \code{\link{combinePValues}}.
#' By default, the contribution from each level is weighted by its number of cells;
#' if \code{equiweight=TRUE}, each level is given equal weight instead.
#' 
#' @author Aaron Lun
#'
#' @seealso
#' \code{\link{fitLinearModel}}, which performs the hard work of fitting the linear models.
#'
#' @examples
#' y <- matrix(rnorm(10000), ncol=100)
#'
#' # Example with categorical factors:
#' A <- gl(2, 50)
#' design <- model.matrix(~A)
#' testLinearModel(y, design, contrasts=c(0, 1))
#'
#' # Example with continuous variables:
#' u <- runif(100)
#' design <- model.matrix(~u)
#' testLinearModel(y, design, contrasts=c(0, 1))
#'
#' # Example with multiple variables:
#' B <- gl(4, 25)
#' design <- model.matrix(~B)
#' testLinearModel(y, design, contrasts=cbind(c(0,1,0,0), c(0,0,1,-1)))
#'
#' @name testLinearModel
NULL

###########################################################

#' @importFrom BiocParallel SerialParam
#' @importFrom beachmat rowBlockApply
#' @importFrom stats p.adjust
.test_linear_model <- function(x, design, coefs=ncol(design), contrasts=NULL, 
    block=NULL, equiweight=FALSE, method="stouffer", subset.row=NULL, BPPARAM=SerialParam()) 
{
    if (!is.null(subset.row)) {
        x <- x[subset.row,,drop=FALSE]
    }

    if (is.null(block)) {
        .test_linear_model_simple(x, design, coefs=coefs, contrasts=contrasts, BPPARAM=BPPARAM)
    } else {
        collated <- split(seq_len(ncol(x)), block)

        # We want the parallelization to be as fine-grained as possible so we
        # do it here; we don't punt it to the fitLinearModel() function, as 
        # then we would have to restart the parallel workers for each block.
        output <- rowBlockApply(x, FUN=.test_linear_model_multiblock, 
            collated=collated, equiweight=equiweight, method=method, 
            design=design, coefs=coefs, contrasts=contrasts, 
            BPPARAM=BPPARAM)

        if (any(vapply(output, is.null, TRUE))) {
            stop("no level of 'block' has a full column rank 'design'")
        }

        output <- do.call(rbind, output)

        # Resetting all the FDRs.
        output$FDR <- p.adjust(output$p.value, method="BH")
        for (i in seq_along(output$per.block)) {
            output$per.block[[i]]$FDR <- p.adjust(output$per.block[[i]]$p.value, method="BH")
        }

        output
    }
}

#' @importFrom S4Vectors metadata
.test_linear_model_multiblock <- function(collated, x, design, equiweight, method, ...) {
    ncells <- lengths(collated)
    for (i in seq_along(collated)) {
        sub <- collated[[i]] 
        res <- .test_linear_model_simple(x[,sub,drop=FALSE], design=design[sub,,drop=FALSE], ...)

        collated[[i]] <- res
        if (is.na(metadata(res)$residual.df) || metadata(res)$residual.df==0L) {
            ncells[i] <- -Inf
        }
    }

    if (all(ncells < 0L)) {
        return(NULL)
    }

    targets <- setdiff(colnames(collated[[1]]), c("p.value", "FDR"))
    output <- combineBlocks(collated, 
        method=method, 
        geometric=FALSE,
        equiweight=equiweight, 
        weights=ncells, 
        ave.fields=targets,
        pval.field="p.value", 
        valid=ncells > 0L)

    rownames(output) <- rownames(collated[[1]])
    output
}

#' @importFrom scuttle fitLinearModel
#' @importFrom limma lmFit classifyTestsF contrasts.fit
#' @importFrom S4Vectors DataFrame metadata<- metadata
#' @importFrom stats p.adjust pt pf
.test_linear_model_simple <- function(x, design, coefs=ncol(design), contrasts=NULL, ...) {
    full <- fitLinearModel(x, design, get.coefs=TRUE, rank.error=FALSE, ...)

    if (is.null(contrasts)) {
        contrasts <- matrix(0, ncol(design), length(coefs))
        contrasts[cbind(coefs, seq_along(coefs))] <- 1
        if (length(coefs) > 1) {
            colnames(contrasts) <- colnames(design)[coefs]
        }
    } else {  
        if (is.null(dim(contrasts))) {
            contrasts <- matrix(contrasts)
        }
    }
    if (ncol(contrasts)==1L && is.null(colnames(contrasts))) {
        colnames(contrasts) <- "logFC"
    }

    if (is.na(full$residual.df)) {
        pvalue <- rep(NA_real_, nrow(full$coefficients))
        coefs <- matrix(NA_real_, length(pvalue), ncol(contrasts))
        colnames(coefs) <- colnames(contrasts)

    } else {
        # Hacking limma to compute our desired statistics.
        lfit <- lmFit(rbind(seq_len(nrow(design))), design)
        lfit$coefficients <- full$coefficients
        lfit$sigma2 <- full$variance
        lfit <- contrasts.fit(lfit, contrasts)

        coefs <- lfit$coefficients
        tstat <- coefs / outer(sqrt(lfit$sigma2), lfit$stdev.unscaled[1,])

        if (ncol(tstat)==1L) {
            tstat <- drop(tstat)
            pvalue <- pt(abs(tstat), df=full$residual.df, lower.tail=FALSE) * 2
        } else {
            lfit$tstat <- tstat
            fstat <- classifyTestsF(lfit, fstat.only=TRUE)
            pvalue <- pf(fstat, ncol(tstat), full$residual.df, lower.tail = FALSE)
            attributes(pvalue) <- NULL
        }
    }

    output <- DataFrame(row.names=rownames(full$coefficients), # account for subsetting.
        coefs, 
        p.value=pvalue,
        FDR=p.adjust(pvalue, method="BH"))

    metadata(output)$residual.df <- full$residual.df
    output
}

###########################################################

#' @export
#' @rdname testLinearModel
setGeneric("testLinearModel", function(x, ...) standardGeneric("testLinearModel"))

#' @export
#' @rdname testLinearModel
setMethod("testLinearModel", "ANY", .test_linear_model)

#' @export
#' @rdname testLinearModel
#' @importFrom SummarizedExperiment assay
setMethod("testLinearModel", "SummarizedExperiment", function(x, ..., assay.type="logcounts") {
    .test_linear_model(assay(x, assay.type), ...)
})
