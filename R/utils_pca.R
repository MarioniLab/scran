#' @importFrom BiocSingular runSVD bsparam
#' @importFrom BiocParallel SerialParam
.centered_SVD <- function(y, max.rank, BSPARAM=bsparam(), BPPARAM=SerialParam(), keep.left=TRUE, keep.right=TRUE)
# Performs the PCA given a log-expression matrix.
# Switches between svd() and irlba() on request.
# Output format is guaranteed to be the same.
{
    runSVD(y, center=TRUE, BSPARAM=BSPARAM, k=max.rank, 
        nu=if (keep.left) max.rank else 0L,
        nv=if (keep.right) max.rank else 0L,
        BPPARAM=BPPARAM)
}

.svd_to_pca <- function(svd.out, ncomp, named=TRUE) 
# Converts centred results to PCs.
{
    if (is.null(svd.out$u)) {
        stop("missing 'U' in SVD results")
    } else if (ncomp > ncol(svd.out$u)) {
        warning("requested number of components greater than available rank")
        ncomp <- ncol(svd.out$u)
    }

    ix <- seq_len(ncomp)
    U <- svd.out$u[,ix,drop=FALSE]
    D <- svd.out$d[ix]

    # Pulling out the PCs (column-multiplying the left eigenvectors).
    pcs <- sweep(U, 2L, D, FUN="*", check.margin = FALSE)
    if (named) {
        colnames(pcs) <- sprintf("PC%i", ix)
    }
    pcs
}

#' @importFrom Matrix rowMeans colSums
.svd_to_rot <- function(svd.out, ncomp, original.mat, subset.row, fill.missing) {
    ncomp <- min(ncomp, ncol(svd.out$v))

    ix <- seq_len(ncomp)
    V <- svd.out$v[,ix,drop=FALSE]
    if (is.null(subset.row) || !fill.missing) {
        rownames(V) <- rownames(original.mat)[subset.row]
        return(V)
    }

    U <- svd.out$u[,ix,drop=FALSE]
    D <- svd.out$d[ix]

    fullV <- matrix(0, nrow(original.mat), ncomp)
    rownames(fullV) <- rownames(original.mat)
    colnames(fullV) <- colnames(V)
    fullV[subset.row,] <- V

    # The idea is that after our SVD, we have X=UDV' where each column of X is a gene.
    # Leftover genes are new columns in X, which are projected on the space of U by doing U'X.
    # This can be treated as new columns in DV', which can be multiplied by U to give denoised values.
    # I've done a lot of implicit transpositions here, hence the code does not tightly follow the logic above.
    leftovers <- !logical(nrow(original.mat))
    leftovers[subset.row] <- FALSE

    left.x <- original.mat[leftovers,,drop=FALSE] 
    left.x <- as.matrix(left.x %*% U) - outer(rowMeans(left.x), colSums(U))

    fullV[leftovers,] <- sweep(left.x, 2, D, "/", check.margin=FALSE)

    fullV
}

#' @importFrom BiocSingular LowRankMatrix
#' @importFrom SingleCellExperiment reduced.dim.matrix reducedDim<-
#' @importFrom SummarizedExperiment assay<-
.pca_to_output <- function(x, pcs, value=c("pca", "lowrank"), name="PCA") { 
    if (value=="pca"){
        out <- reduced.dim.matrix(pcs$components)
        attr(out, "percentVar") <- pcs$percent.var
        attr(out, "varExplained") <- pcs$var.explained
        attr(out, "rotation") <- pcs$rotation
    } else {
        out <- LowRankMatrix(pcs$rotation, pcs$components)
    }

    value <- match.arg(value)
    if (value=="pca"){
        if (is.null(name)) name <- "PCA"
        reducedDim(x, name) <- out
    } else if (value=="lowrank") {
        if (is.null(name)) name <- "lowrank"
        assay(x, i=name) <- out
    }
    x
}

#' @importFrom scuttle .subset2index
.process_subset_for_pca <- function(subset.row, x) {
    if (missing(subset.row)) {
        warning(paste(strwrap("'subset.row=' is typically used to specify HVGs for PCA. If the use of all genes is intentional, suppress this message with 'subset.row=NULL'."), collapse="\n"))
        subset.row <- NULL
    }
    .subset2index(subset.row, x, byrow=TRUE)
}
