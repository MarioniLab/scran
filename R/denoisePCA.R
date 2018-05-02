#' @importFrom DelayedArray DelayedArray
#' @importFrom DelayedMatrixStats rowVars rowMeans2
.denoisePCA <- function(x, technical, subset.row=NULL,
                        value=c("pca", "n", "lowrank"), min.rank=5, max.rank=100, 
                        approximate=FALSE, rand.seed=1000, irlba.args=list())
# Performs PCA and chooses the number of PCs to keep based on the technical noise.
# This is done on the residuals if a design matrix is supplied.
#
# written by Aaron Lun
# created 13 March 2017    
{
    subset.row <- .subset_to_index(subset.row, x, byrow=TRUE)
    x2 <- DelayedArray(x)
    scale <- NULL

    # Processing different mechanisms through which we specify the technical component.
    if (is(technical, "DataFrame")) { 
        x.var <- rowVars(x2)
        scale <- sqrt(technical$total/x.var) # Making sure everyone has the reported total variance.
        all.var <- technical$total[subset.row]
        tech.var <- technical$tech[subset.row]
    } else {
        all.var <- rowVars(x2, rows=subset.row)
        if (is.function(technical)) {
            all.means <- rowMeans2(x2, rows=subset.row)
            tech.var <- technical(all.means)
        } else {
            tech.var <- technical[subset.row]
        }
    }

    # Filtering out genes with negative biological components.
    keep <- all.var > tech.var
    use.rows <- subset.row[keep]
    tech.var <- tech.var[keep]
    all.var <- all.var[keep]
    total.tech <- sum(tech.var)

    # Subsetting and scaling the matrix.
    y <- x[use.rows,,drop=FALSE] 
    if (!is.null(scale)) {
        y <- y * scale[use.rows]
    }

    # Setting up the PCA function and its arguments.
    value <- match.arg(value)
    args <- list(y=t(y), max.rank=max.rank, value=value)
    if (approximate) {
        svdfun <- .irlba_svd
        args <- c(args, irlba.args)
    } else {
        svdfun <- .full_svd
    }

    # Runing the PCA and choosing the number of PCs.
    original <- do.call(svdfun, args)
    var.exp <- original$d^2 / (ncol(y) - 1)
    npcs <- .get_npcs_to_keep(var.exp, total.tech, total=sum(all.var))
    npcs <- .keep_rank_in_range(npcs, min.rank, length(var.exp))

    # Processing remaining aspects.
    out.val <- .convert_to_output(original, npcs, value, x, scale, use.rows)
    attr(out.val, "percentVar") <- var.exp/sum(all.var)
    return(out.val)
} 

.get_npcs_to_keep <- function(var.exp, technical, total=sum(var.exp)) 
# Discarding PCs until we get rid of as much technical noise as possible
# while preserving the biological signal. This is done by assuming that 
# the biological signal is fully contained in earlier PCs, such that we 
# discard the later PCs until we account for 'tech.var'.
{
    npcs <- length(var.exp)
    flipped.var.exp <- rev(var.exp)
    estimated.contrib <- cumsum(flipped.var.exp) 

    above.noise <- estimated.contrib > technical - (total - sum(var.exp)) 
    if (any(above.noise)) { 
        to.keep <- npcs - min(which(above.noise)) + 1L
    } else {
        to.keep <- npcs
    }
    return(to.keep)
}

##############################
# S4 method definitions here #
##############################

#' @export
setGeneric("denoisePCA", function(x, ...) standardGeneric("denoisePCA"))

#' @export
setMethod("denoisePCA", "ANY", .denoisePCA)

#' @importFrom SummarizedExperiment assay "assay<-"
#' @importFrom SingleCellExperiment reducedDim isSpike
#' @export
setMethod("denoisePCA", "SingleCellExperiment", 
          function(x, ..., subset.row=NULL, value=c("pca", "n", "lowrank"), 
                   assay.type="logcounts", get.spikes=FALSE, sce.out=TRUE) {

    subset.row <- .SCE_subset_genes(subset.row=subset.row, x=x, get.spikes=get.spikes)
    out <- .denoisePCA(assay(x, i=assay.type), ..., value=value, subset.row=subset.row)

    value <- match.arg(value) 
    if (!sce.out || value=="n") { 
        return(out)
    }

    if (value=="pca"){ 
        reducedDim(x, "PCA") <- out
    } else if (value=="lowrank") {
        if (!get.spikes) {
            out[isSpike(x),] <- 0
        }
        assay(x, i="lowrank") <- out
    }
    return(x)
})

