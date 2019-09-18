#' Normalization by deconvolution
#'
#' Scaling normalization of single-cell RNA-seq data by deconvolving size factors from cell pools.
#' 
#' @param x For \code{calculateSumFactors}, a  numeric matrix-like object of counts, where rows are genes and columns are cells.
#' Alternatively, a \linkS4class{SummarizedExperiment} object containing such a matrix.
#'
#' For \code{computeSumFactors}, a \linkS4class{SingleCellExperiment} object containing a count matrix.
#' @param sizes A numeric vector of pool sizes, i.e., number of cells per pool.
#' @param clusters An optional factor specifying which cells belong to which cluster, for deconvolution within clusters.
#' @param ref.clust A level of \code{clusters} to be used as the reference cluster for inter-cluster normalization.
#' @param max.cluster.size An integer scalar specifying the maximum number of cells in each cluster.
#' @param positive A logical scalar indicating whether linear inverse models should be used to enforce positive estimates.
#' @param scaling A numeric scalar containing scaling factors to adjust the counts prior to computing size factors.
#' @param min.mean A numeric scalar specifying the minimum (library size-adjusted) average count of genes to be used for normalization.
#' @param subset.row See \code{?"\link{scran-gene-selection}"}.
#' @param BPPARAM A BiocParallelParam object specifying whether and how clusters should be processed in parallel.
#' @param ... For the \code{calculateSumFactors} generic, additional arguments to pass to each method.
#' For the \linkS4class{SummarizedExperiment} method, additional methods to pass to the ANY method.
#' 
#' For the \code{computeSumFactors} function, additional arguments to pass to \code{calculateSumFactors}.
#' @param assay.type A string specifying which assay values to use when \code{x} is a SummarizedExperiment or SingleCellExperiment.
#' @param get.spikes See \code{?"\link{scran-gene-selection}"}.
#' @param sf.out Deprecated, a logical scalar indicating whether only size factors should be returned.
#' 
#' @section Overview of the deconvolution method:
#' The \code{computeSumFactors} function implements the deconvolution strategy (Lun et al., 2016) for scaling normalization of sparse count data.
#' Briefly, a pool of cells is selected and the expression profiles for those cells are summed together.
#' The pooled expression profile is normalized against an average reference pseudo-cell, constructed by averaging the counts across all cells.
#' This defines a size factor for the pool as the median ratio between the count sums and the average across all genes.
#' 
#' The scaling bias for the pool is equal to the sum of the biases for the constituent cells.
#' The same applies for the size factors, as these are effectively estimates of the bias for each cell.
#' This means that the size factor for the pool can be written as a linear equation of the size factors for the cells.
#' Repeating this process for multiple pools will yield a linear system that can be solved to obtain the size factors for the individual cells.
#' 
#' In this manner, pool-based factors are deconvolved to yield the relevant cell-based factors.
#' The advantage is that the pool-based estimates are more accurate, as summation reduces the number of stochastic zeroes and the associated bias of the size factor estimate.
#' This accuracy feeds  back into the deconvolution process, thus improving the accuracy of the cell-based size factors.
#' 
#' @section Pooling with a sliding window:
#' Within each cluster (if not specified, all cells are put into a single cluster), cells are sorted by increasing library size and a sliding window is applied to this ordering.
#' Each location of the window defines a pool of cells with similar library sizes.
#' This avoids inflated estimation errors for very small cells when they are pooled with very large cells.
#' Sliding the window will construct an over-determined linear system that can be solved by least-squares methods to obtain cell-specific size factors.
#' 
#' Window sliding is repeated with different window sizes to construct the linear system, as specified by \code{sizes}.
#' By default, the number of cells in each window ranges from 21 to 101.
#' Using a range of window sizes improves the precision of the estimates, at the cost of increased computational complexity.
#' The defaults were chosen to provide a reasonable compromise between these two considerations.
#' The default choice also avoids rare cases of linear dependencies and unstable estimates when all pool sizes are not co-prime with the number of cells.
#' 
#' The smallest window should be large enough so that the pool-based size factors are, on average, non-zero.
#' We recommend window sizes no lower than 20 for UMI data, though smaller windows may be possible for read count data.
#' 
#' If there are fewer cells than the smallest window size, the function will naturally degrade to performing library size normalization.
#' This yields results that are the same as \code{\link{librarySizeFactors}}.
#' 
#' @section Prescaling to improve accuracy:
#' The simplest approach to pooling is to simply add the counts together for all cells in each pool.
#' However, this is suboptimal as any errors in the estimation of the pooled size factor will propagate to all component cell-specific size factors upon solving the linear system.
#' If the error is distributed evenly across all cell-specific size factors, the small size factors will have larger relative errors compared to the large size factors.
#' 
#' To avoid this, we perform \dQuote{prescaling} where we divide the counts by a cell-specific factor prior to pooling.
#' The prescaling factor should be close to the true size factor for each cell.
#' (Obviously, the true size factor is unknown so we use the library size for each cell instead.)
#' Solving the linear system constructed with prescaled values should yield estimates that are more-or-less equal across all cells.
#' Thus, given similar absolute errors, the relative errors for all cells will also be similar.
#' 
#' The default usage assumes that the library size is roughly proportional to the true size factor.
#' This can be violated in pathological scenarios involving extreme differential expression.
#' In such cases, we recommend running \code{computeSumFactors} twice to improve accuracy.
#' The first run is done as usual and will yield an initial estimate of the size factor for each cell.
#' In the second run, we use our initial estimates for prescaling by supplying them to the \code{scaling} argument, 
#' This weakens the assumption above by avoiding large inaccuracies in the prescaling factor.
#' 
#' Obviously, this involves twice as much computational work.
#' As a result, performing multiple iterations is not the default recommendation, especially as the benefits are only apparent in extreme circumstances.
#' 
#' @section Solving the linear system:
#' The linear system is solved using the sparse QR decomposition from the \pkg{Matrix} package.
#' However, this has known problems when the linear system becomes too large (see \url{https://stat.ethz.ch/pipermail/r-help/2011-August/285855.html}).
#' In such cases, we set \code{clusters} to break up the linear system into smaller, more manageable components that can be solved separately.
#' The default \code{max.cluster.size} will arbitrarily break up the cell population (within each cluster, if specified) so that we never pool more than 3000 cells.
#' 
#' @section Normalization within and between clusters:
#' In general, it is more appropriate to pool more similar cells to avoid violating the assumption of a non-DE majority of genes across the data set.
#' This can be done by specifying the \code{clusters} argument where cells in each cluster have similar expression profiles.
#' Deconvolution is subsequently applied on the cells within each cluster.
#' A convenience function \code{\link{quickCluster}} is provided for this purpose, though any reasonable clustering can be used.
#' 
#' Size factors computed within each cluster must be rescaled for comparison between clusters.
#' This is done by normalizing between clusters to identify the rescaling factor.
#' One cluster is chosen as a ``reference'' to which all others are normalized.
#' Ideally, the reference cluster should have a stable expression profile and not be extremely different from all other clusters.
#' 
#' By default, the cluster with the most non-zero counts is used as the reference.
#' This reduces the risk of obtaining undefined rescaling factors for the other clusters, while improving the precision (and also accuracy) of the median-based estimate of each factor.
#' Alternatively, the reference can be manually specified using \code{ref.clust} if there is prior knowledge about which cluster is most suitable, e.g., from PCA or t-SNE plots.
#' 
#' Each cluster should ideally be large enough to contain a sufficient number of cells for pooling.
#' Otherwise, \code{computeSumFactors} will degrade to library size normalization.
#' 
#' @section Dealing with negative size factors:
#' In theory, it is possible to obtain negative estimates for the size factors.
#' These values are obviously nonsensical and \code{computeSumFactors} will raise a warning if they are encountered.
#' Negative estimates are mostly commonly generated from low quality cells with few expressed features, such that most counts are zero even after pooling.
#' They may also occur if insufficient filtering of low-abundance genes was performed.
#' 
#' To avoid negative size factors, the best solution is to increase the stringency of the filtering.
#' \itemize{
#' \item If only a few negative size factors are present, they are likely to correspond to a few low-quality cells with few expressed features.
#' Such cells are difficult to normalize reliably under any approach, and can be removed by increasing the stringency of the quality control.
#' \item If many negative size factors are present, it is probably due to insufficient filtering of low-abundance genes.
#' This results in many zero counts and pooled size factors of zero, and can be fixed by filtering out more genes with a higher \code{min.mean}.
#' }
#' Another approach is to increase in the number of \code{sizes} to improve the precision of the estimates.
#' This reduces the chance of obtaining negative size factors due to estimation error, for cells where the true size factors are very small.
#' 
#' As a last resort, \code{positive=TRUE} is set by default, which uses \code{\link{cleanSizeFactors}} to coerce any negative estimates to positive values.
#' This ensures that, at the very least, downstream analysis is possible even if the size factors for affected cells are not accurate.
#' Users can skip this step by setting \code{positive=FALSE} to perform their own diagnostics or coercions.
#' 
#' @section Gene selection:
#' If too many genes have consistently low counts across all cells, even the pool-based size factors will be zero.
#' This results in zero or negative size factor estimates for many cells.
#' We avoid this by filtering out low-abundance genes using the \code{min.mean} argument.
#' This represents a minimum threshold \code{min.mean} on the library size-adjusted average counts from \code{\link{calcAverage}}.
#' 
#' By default, we set \code{min.mean} to 1 for read count data and 0.1 for UMI data.
#' The exact values of these defaults are more-or-less arbitrary and are retained for historical reasons.
#' The lower threshold for UMIs is motivated by (i) their lower count sizes, which would result in the removal of too many genes with a higher threshold; and (ii) the lower variability of UMI counts, which results in a lower frequency of zeroes compared to read count data at the same mean.
#' We use the mean library size to detect whether the counts are those of reads (above 100,000) or UMIs (below 50,000) to automatically set \code{min.mean}.
#' Mean library sizes in between these two limits will trigger a warning and revert to using \code{min.mean=0.1}.
#' 
#' If \code{clusters} is specified, filtering by \code{min.mean} is performed on the per-cluster average during within-cluster normalization,
#' and then on the (library size-adjusted) average of the per-cluster averages during between-cluster normalization.
#' 
#' @section Obtaining standard errors:
#' Previous versions of \code{computeSumFactors} would return the standard error for each size factor when \code{errors=TRUE}.
#' This argument is no longer available as we have realized that standard error estimation from the linear model is not reliable.
#' Errors are likely underestimated due to correlations between pool-based size factors when they are computed from a shared set of underlying counts.
#' Users wishing to obtain a measure of uncertainty are advised to perform simulations instead, using the original size factor estimates to scale the mean counts for each cell.
#' Standard errors can then be calculated as the standard deviation of the size factor estimates across simulation iterations.
#' 
#' @return
#' For \code{calculateSumFactors}, a numeric vector of size factors for all cells in \code{x} is returned.
#' 
#' For \code{computeSumFactors}, an object of class \code{x} is returned containing the vector of size factors in \code{\link{sizeFactors}(x)}.
#' 
#' @author
#' Aaron Lun and Karsten Bach
#' 
#' @seealso
#' \code{\link{quickCluster}}, to obtain a rough clustering for use in \code{clusters}.
#'
#' \code{\link{logNormCounts}}, which uses the computed size factors to compute normalized expression values.
#'
#' @examples
#' library(scater)
#' sce <- mockSCE(ncells=500)
#' 
#' # Computing the size factors.
#' sce <- computeSumFactors(sce)
#' head(sizeFactors(sce))
#' plot(librarySizeFactors(sce), sizeFactors(sce), log="xy")
#'
#' # Using pre-clustering.
#' preclusters <- quickCluster(sce)
#' table(preclusters)
#' 
#' sce2 <- computeSumFactors(sce, clusters=preclusters)
#' head(sizeFactors(sce2))
#'
#' @references
#' Lun ATL, Bach K and Marioni JC (2016).
#' Pooling across cells to normalize single-cell RNA sequencing data with many zero counts.
#' \emph{Genome Biol.} 17:75
#'
#' @name computeSumFactors
NULL

#' @importFrom BiocParallel bplapply SerialParam
.calculate_sum_factors <- function(x, sizes=seq(21, 101, 5), clusters=NULL, ref.clust=NULL, max.cluster.size=3000, 
    positive=TRUE, scaling=NULL, min.mean=NULL, subset.row=NULL, BPPARAM=SerialParam())
# This contains the function that performs normalization on the summed counts.
# It also provides support for normalization within clusters, and then between
# clusters to make things comparable. 
{
    ncells <- ncol(x)
    if (is.null(clusters)) {
        clusters <- integer(ncells)
    }
	clusters <- .limit_cluster_size(clusters, max.cluster.size)

    if (ncells!=length(clusters)) { 
        stop("'ncol(x)' is not equal to 'length(clusters)'")
    }
    indices <- split(seq_along(clusters), clusters)

    if (length(indices)==0L || any(lengths(indices)==0L)) {
        stop("zero cells in one of the clusters")
    }

    # Checking sizes and subsetting.
    sizes <- sort(as.integer(sizes))
    if (anyDuplicated(sizes)) { 
        stop("'sizes' are not unique") 
    }
    subset.row <- .subset_to_index(subset.row, x, byrow=TRUE)

    # Setting some other values.
    nclusters <- length(indices)
    clust.nf <- clust.profile <- clust.libsizes <- vector("list", nclusters)
    clust.meanlib <- numeric(nclusters)

    # Computing normalization factors within each cluster.
    all.norm <- bplapply(indices, FUN=.per_cluster_normalize, x=x, sizes=sizes, subset.row=subset.row, 
        min.mean=min.mean, positive=positive, scaling=scaling, BPPARAM=BPPARAM)

    clust.nf <- lapply(all.norm, "[[", i="final.nf")
    clust.profile <- lapply(all.norm, "[[", i="ave.cell")

    # Adjusting size factors between clusters.
    if (is.null(ref.clust)) {
        non.zeroes <- vapply(clust.profile, FUN=function(x) sum(x>0), FUN.VALUE=0L) 
        ref.clust <- which.max(non.zeroes)
    }
    rescaling.factors <- .rescale_clusters(clust.profile, ref.col=ref.clust, min.mean=min.mean) 

    clust.nf.scaled <- vector("list", nclusters)
    for (clust in seq_len(nclusters)) { 
        clust.nf.scaled[[clust]] <- clust.nf[[clust]] * rescaling.factors[[clust]]
    }
    clust.nf.scaled <- unlist(clust.nf.scaled)

    # Returning centered size factors, rather than normalization factors.
    final.sf <- rep(NA_integer_, ncells)
    indices <- unlist(indices)
    final.sf[indices] <- clust.nf.scaled
    
    is.pos <- final.sf > 0 & !is.na(final.sf)
    final.sf/mean(final.sf[is.pos])
}

#############################################################
# Internal functions.
#############################################################

#' @importFrom Matrix qr qr.coef colSums
#' @importFrom scater nexprs
.per_cluster_normalize <- function(x, curdex, sizes, subset.row, min.mean=NULL, positive=FALSE, scaling=NULL) 
# Computes the normalization factors _within_ each cluster,
# along with the reference pseudo-cell used for normalization. 
# Written as a separate function so that bplapply operates in the scran namespace.
{
    cur.cells <- length(curdex)
    sizes <- sizes[sizes <= cur.cells]

    vals <- subset_and_divide(x, subset.row-1L, curdex-1L, scaling)
    scaling <- vals[[1]]
    exprs <- vals[[2]]
    ave.cell <- vals[[3]] # equivalent to calcAverage().

    # Choosing a mean filter based on the data type and then filtering:
    if (is.null(min.mean)) {
        mean.lib <- sum(ave.cell)
        if (mean.lib <= 50000) { # Probably UMI data.
            min.mean <- 0.1
        } else if (mean.lib >= 100000) { # Probably read data.
            min.mean <- 1
        } else {
            min.mean <- 0.1
            warning("assuming UMI data when setting 'min.mean'")
        }
    }

    min.mean <- pmax(min.mean, 1e-8) # must be positive.
    high.ave <- min.mean <= ave.cell 
    use.ave.cell <- ave.cell
    if (!all(high.ave)) { 
        exprs <- exprs[high.ave,,drop=FALSE]
        use.ave.cell <- use.ave.cell[high.ave]
    }

    # Using our summation approach.
    sphere <- .generateSphere(scaling)
    new.sys <- .create_linear_system(exprs, use.ave.cell, sphere, sizes) 
    design <- new.sys$design
    output <- new.sys$output

    # Weighted least-squares.
    QR <- qr(design)
    final.nf <- qr.coef(QR, output)
    final.nf <- final.nf * scaling

    if (any(final.nf < 0)) {
        warning("encountered negative size factor estimates")
        if (positive) {
            num.detected <- nexprs(exprs, byrow=FALSE)
            final.nf <- cleanSizeFactors(final.nf, num.detected) 
        }
    }

    list(final.nf=final.nf, ave.cell=ave.cell)
}

.generateSphere <- function(lib.sizes) 
# Sorts cells by their library sizes, and generates an ordering vector
# to arrange cells in a circle based on increasing/decreasing lib size.
{
    nlibs <- length(lib.sizes)
    o <- order(lib.sizes)
    even <- seq(2,nlibs,2)
    odd <- seq(1,nlibs,2)
    out <- c(o[odd], rev(o[even]))
    c(out, out)
}

LOWWEIGHT <- 0.000001

#' @importFrom Matrix sparseMatrix
.create_linear_system <- function(cur.exprs, ave.cell, sphere, pool.sizes) 
# Does the heavy lifting of computing pool-based size factors 
# and creating the linear system out of the equations for each pool.
{
    row.dex <- col.dex <- output <- vector("list", 2L)

    # Creating the linear system with the requested pool sizes.
    out <- pool_size_factors(cur.exprs, ave.cell, sphere - 1L, pool.sizes)
    row.dex[[1]] <- out[[1]] + 1L
    col.dex[[1]] <- out[[2]] + 1L
    output[[1]]<- out[[3]]

    # Adding extra equations to guarantee solvability.
    cur.cells <- ncol(cur.exprs)
    row.dex[[2]] <- seq_len(cur.cells) + cur.cells * length(pool.sizes)
    col.dex[[2]] <- seq_len(cur.cells)
    output[[2]] <- rep(sqrt(LOWWEIGHT) / sum(ave.cell), cur.cells) # equivalent to library size factors for each cell, but downweighted.

    # Setting up the entries of the LHS matrix.
    eqn.values <- rep(c(1, sqrt(LOWWEIGHT)), lengths(row.dex))

    # Constructing a sparse matrix.
    row.dex <- unlist(row.dex)
    col.dex <- unlist(col.dex)
    output <- unlist(output)
    design <- sparseMatrix(i=row.dex, j=col.dex, x=eqn.values, dims=c(length(output), cur.cells))

    return(list(design=design, output=output))
}

#' @importFrom stats median
.rescale_clusters <- function(mean.prof, ref.col, min.mean) 
# Chooses a cluster as a reference and rescales all other clusters to the reference,
# based on the 'normalization factors' computed between pseudo-cells.
{
    if (is.character(ref.col)) {
        ref.col <- which(names(mean.prof)==ref.col)
        if (length(ref.col)==0L) { 
            stop("'ref.clust' not in 'clusters'")
        }
    }

    nclusters <- length(mean.prof)
    rescaling <- numeric(nclusters)
    for (clust in seq_len(nclusters)) { 
        ref.prof <- mean.prof[[ref.col]]
        cur.prof <- mean.prof[[clust]] 

        # Filtering based on the mean of the per-cluster means (requires scaling for the library size).
        # Effectively equivalent to 'calcAverage(cbind(ref.ave.count, cur.ave.count))' where the averages
        # are themselves equivalent to 'calcAverage()' across all cells in each cluster.
        cur.libsize <- sum(cur.prof)
        ref.libsize <- sum(ref.prof)
        to.use <- (cur.prof/cur.libsize + ref.prof/ref.libsize)/2 * (cur.libsize + ref.libsize)/2 >= min.mean
        if (!all(to.use)) { 
            cur.prof <- cur.prof[to.use]
            ref.prof <- ref.prof[to.use]
        } 

        # Adjusting for systematic differences between clusters.
        rescaling[[clust]] <- median(cur.prof/ref.prof, na.rm=TRUE)
    }

    if (any(!is.finite(rescaling) | rescaling<=0)) {
        stop("inter-cluster rescaling factors are not strictly positive")
    }
    names(rescaling) <- names(mean.prof)
    rescaling
}

.limit_cluster_size <- function(clusters, max.size) 
# Limits the maximum cluster size to avoid problems with memory in Matrix::qr().
# Done by arbitrarily splitting large clusters so that they fall below max.size.
{
    if (is.null(max.size)) { 
        return(clusters) 
    }
    
    new.clusters <- integer(length(clusters))
    counter <- 1L
    for (id in unique(clusters)) {
        current <- id==clusters
        ncells <- sum(current)
        
        if (ncells <= max.size) {
            new.clusters[current] <- counter
            counter <- counter+1L
            next
        }
       
        # Size of output clusters is max.size * N / ceil(N), where N = ncells/max.size.
        # This is minimal at the smallest N > 1, where output clusters are at least max.size/2. 
        # Thus, we need max.size/2 >= min.size to guarantee that the output clusters are >= min.size.
        mult <- ceiling(ncells/max.size)
        realloc <- rep(seq_len(mult) - 1L + counter, length.out=ncells)
        new.clusters[current] <- realloc
        counter <- counter + mult
    }

    factor(new.clusters)
}

#############################################################
# S4 method definitions.
#############################################################

#' @export
#' @rdname computeSumFactors
setGeneric("calculateSumFactors", function(x, ...) standardGeneric("calculateSumFactors"))

#' @export
#' @rdname computeSumFactors
setMethod("calculateSumFactors", "ANY", .calculate_sum_factors)

#' @export
#' @rdname computeSumFactors
#' @importFrom SummarizedExperiment assay
setMethod("calculateSumFactors", "SummarizedExperiment", function(x, ..., assay.type="counts") {
    .calculate_sum_factors(assay(x, i=assay.type), ...)
})

#' @export
#' @rdname computeSumFactors
#' @importFrom SummarizedExperiment assay 
#' @importFrom BiocGenerics "sizeFactors<-"
computeSumFactors <- function(x, ..., subset.row=NULL, assay.type="counts", get.spikes=FALSE, sf.out=FALSE) { 
    if (!is(x, "SingleCellExperiment")) {
        .Deprecated(msg="use 'calculateSumFactors' for any 'x' that is not a SingleCellExperiment")
        return(.calculate_sum_factors(x, ..., subset.row=subset.row))
    }

    subset.row <- .SCE_subset_genes(subset.row=subset.row, x=x, get.spikes=get.spikes)
    sf <- .calculate_sum_factors(assay(x, i=assay.type), subset.row=subset.row, ...) 
    if (sf.out) { 
        .Deprecated(old="'sf.out=TRUE'", new="calculateSumFactors")
        return(sf) 
    }
    sizeFactors(x) <- sf
    x
}
