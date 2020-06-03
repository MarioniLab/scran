#' Test for differences along pseudotime
#'
#' Implements a simple method of testing for significant differences with respect to pseudotime,
#' based on fitting linear models with a spline basis matrix.
#'
#' @param x A numeric matrix-like object containing log-expression values for cells (columns) and genes (rows).
#' Alternatively, a \linkS4class{SummarizedExperiment} containing such a matrix.
#' @param pseudotime A numeric matrix with one row per cell in \code{x} and one column per path (i.e., lineage).
#' A vector is treated the same as a 1-column matrix.
#' @param df Integer scalar specifying the degrees of freedom for the splines.
#' @param get.lfc Logical scalar indicating whether to return an overall log-fold change along each path.
#' @param get.spline.coef Logical scalar indicating whether to return the estimates of the spline coefficients.
#' @param trend.only Logical scalar indicating whether only differences in the trend should be considered
#' when testing for differences between paths.
#' @param ... For the generic, further arguments to pass to specific method.
#'
#' For the ANY method, further arguments to pass to \code{\link{fitLinearModel}}.
#' 
#' For the SummarizedExperiment method, further arguments to pass to the ANY method.
#' @param assay.type String or integer scalar specifying the assay containing the log-expression matrix.
#'
#' @return
#' A \linkS4class{DataFrame} is returned containing the statistics for each gene (row),
#' including the p-value and its BH-adjusted equivalent.
#' If \code{get.lfc=TRUE}, an overall log-fold change is returned for each path.
#'
#' If \code{get.spline.coef=TRUE}, the estimated spline coefficients are also returned (single path)
#' or the differences in the spline fits to the first path are returned (multiple paths).
#' 
#' @details
#' For a single path in \code{pseudotime},
#' this function fits a natural spline to the expression of each gene with respect to pseudotime.
#' It then does an ANOVA to test whether any of the spline coefficients are non-zero.
#' In this manner, genes exhibiting a significant (and potentially non-linear) trend
#' with respect to the pseudotime can be detected as those with low p-values.
#' 
#' For multiple paths in \code{pseudotime}, 
#' the null hypothesis is that all paths have the same trend (if \code{trend.only=TRUE})
#' or the same trend and intercept (if \code{FALSE}).
#' This is done by effectively fitting a separate trend to each path 
#' and performing an ANOVA to detect differences in the trend alone or in the trend and intercept.
#' In this manner, genes exhibiting differences in behavior between paths can be detected.
#'
#' The expected format of \code{pseudotime} is the same as that returned by \code{\link{orderClusterMST}}.
#' Each cell is assigned to a path if it has a non-\code{NA} value in the corresponding column.
#' For single path testing, cells with \code{NA} values in \code{pseudotime} are ignored;
#' for multiple path testing, cells assigned to multiple paths are ignored.
#'
#' By default, estimates of the spline coefficients are not returned as they are difficult to interpret.
#' Rather, a log-fold change of expression along each path is estimated
#' to provide some indication of the overall magnitude and direction of any change.
#' 
#' @author Aaron Lun
#'
#' @examples
#' y <- matrix(rnorm(10000), ncol=100)
#'
#' # Testing for a difference in a single path:
#' u <- runif(100)
#' testPseudotime(y, u)
#'
#' # Testing for differences in multiple paths
#' # by mocking up a pseudotime matrix.
#' p <- cbind(path1=u, path2=u)
#' path1 <- rbinom(length(u), 1, 0.5)==0
#' p[!path1,1] <- NA
#' p[path1,2] <- NA
#' 
#' testPseudotime(y, p)
#'
#' @seealso
#' \code{\link{orderClusterMST}}, to generate the pseudotime matrix.
#'
#' \code{\link{testLinearModel}}, which performs the tests under the hood.
#'
#' @name testPseudotime
NULL

.test_pseudotime <- function(x, pseudotime, df=5, get.lfc=TRUE, get.spline.coef=FALSE, trend.only=TRUE) {
    if (is.null(dim(pseudotime)) || ncol(pseudotime)==1) {
        .test_solo_pseudotime(x, pseudotime, df=df, get.lfc=get.lfc, 
            get.spline.coef=get.spline.coef)
    } else {
        .test_multi_pseudotime(x, pseudotime, df=df, get.lfc=get.lfc, 
            get.spline.coef=get.spline.coef, trend.only=trend.only)
    }
}

#' @importFrom stats model.matrix
#' @importFrom scuttle fitLinearModel
.test_solo_pseudotime <- function(x, pseudotime, df, get.lfc, get.spline.coef) {
    pseudotime <- drop(pseudotime)

    keep <- !is.na(pseudotime)
    pseudotime <- pseudotime[keep]
    x <- x[,keep,drop=FALSE] 
    design <- .forge_spline_basis_design(pseudotime, df=df)
    output <- testLinearModel(x, design=design, coefs=2:ncol(design))
 
    if (get.lfc) {
        prior <- colnames(output)
        design.lfc <- model.matrix(~pseudotime)
        output$logFC <- fitLinearModel(x, design=design.lfc)$coefficients[,2]
        output <- output[,c("logFC", prior)] 
    }

    if (!get.spline.coef) {
        output <- output[,setdiff(colnames(output), colnames(design))]
    }

    output
}

#' @importFrom stats model.matrix
.forge_spline_basis_design <- function(p, df) { 
    # Uniquify'ing to avoid non-full rank problems when
    # many of the quantiles are stacked at the same position. 
    up <- unique(p)
    if (length(up) <= df) {
        stop("'not enough unique pseudotime values for the specified 'df'")
    }

    basis <- splines::ns(up, df=df)
    colnames(basis) <- sprintf("spline%i", seq_len(df))
    basis <- basis[match(p, up),,drop=FALSE]

    cbind(Intercept=rep(1, length(p)), basis)
}

#' @importFrom stats model.matrix
#' @importFrom scuttle fitLinearModel
.test_multi_pseudotime <- function(x, pseudotime, df, get.lfc, get.spline.coef, trend.only) {
    if (anyDuplicated(colnames(pseudotime))) {
        warning("'pseudotime' has duplicated column names")
        colnames(pseudotime) <- NULL
    }
    if (is.null(colnames(pseudotime))) {
        colnames(pseudotime) <- sprintf("path%i", seq_len(ncol(pseudotime)))
    }

    nonna <- !is.na(pseudotime)
    solo <- rowSums(nonna) == 1
    if (!all(solo)) {
        x <- x[,solo,drop=FALSE]
        pseudotime <- pseudotime[solo,,drop=FALSE]
    }

    common <- rowMeans(pseudotime, na.rm=TRUE)
    design.raw <- .forge_spline_basis_design(common, df=df)

    all.x <- all.design <- all.lfc <- list()
    for (i in colnames(pseudotime)) {
        current.pseudo <- pseudotime[,i]
        keep <- !is.na(current.pseudo) 
        current.pseudo <- current.pseudo[keep]

        all.x[[i]] <- x[,keep,drop=FALSE] 
        cur.design <- design.raw[keep,,drop=FALSE] 
        colnames(cur.design) <- paste0(colnames(cur.design), ".", i)
        all.design[[i]] <- cur.design

        if (get.lfc) {
            design.lfc <- model.matrix(~current.pseudo)
            all.lfc[[paste0("logFC.", i)]] <- fitLinearModel(all.x[[i]], design=design.lfc)$coefficients[,2]
        }
    }

    # Creating the design matrix.
    x2 <- do.call(cbind, all.x)
    for (i in seq_along(all.design)) {
        copy <- rep(all.design[i], length(all.design))
        for (j in setdiff(seq_along(copy), i)) {
            copy[[j]][] <- 0
        }
        all.design[[i]] <- do.call(cbind, copy)
    }
    design <- do.call(rbind, all.design)

    # Building up a mega contrast of doom.
    contrastor <- list()
    for (i in seq_len(ncol(pseudotime) - 1)) {
        N <- df + 1L

        con <- matrix(0, ncol(design), N)
        diag(con) <- -1
        con[cbind(i * N + seq_len(N), seq_len(N))] <- 1
        if (trend.only) {
            con <- con[,-1,drop=FALSE]
        }

        colnames(con) <- sprintf("spline%i.%s", seq_len(ncol(con)), colnames(pseudotime)[i+1])
        contrastor[[i]] <- con
    }

    contrast <- do.call(cbind, contrastor)
    output <- testLinearModel(x2, design=design, contrasts=contrast)

    # Organizing the output.
    if (!get.spline.coef) {
        output <- output[,setdiff(colnames(output), colnames(contrast))]
    }
    if (get.lfc) {
        output <- cbind(DataFrame(all.lfc, row.names=rownames(output)), output)
    }

    output
}

#' @export
#' @rdname testPseudotime
setGeneric("testPseudotime", function(x, ...) standardGeneric("testPseudotime"))

#' @export
#' @rdname testPseudotime
setMethod("testPseudotime", "ANY", .test_pseudotime)

#' @export
#' @rdname testPseudotime
#' @importFrom SummarizedExperiment assay
setMethod("testPseudotime", "SummarizedExperiment", function(x, ..., assay.type="logcounts") {
    .test_pseudotime(assay(x, assay.type), ...)
})
