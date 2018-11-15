#' @importFrom BiocParallel bplapply SerialParam
.computeSumFactors <- function(x, sizes=seq(21, 101, 5), clusters=NULL, ref.clust=NULL, max.cluster.size=3000, 
    positive=TRUE, scaling=NULL, min.mean=1, subset.row=NULL, BPPARAM=SerialParam())
# This contains the function that performs normalization on the summed counts.
# It also provides support for normalization within clusters, and then between
# clusters to make things comparable. It can also switch to linear inverse models
# to ensure that the estimates are non-negative.
#
# written by Aaron Lun
# created 23 November 2015
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

    if (is.null(min.mean)) {
        stop("set 'min.mean=0' to turn off abundance filtering")
    }
    min.mean <- pmax(min.mean, 1e-8) # must be at least non-zero mean. 

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
    final.sf <- final.sf/mean(final.sf[is.pos])
    return(final.sf)
}

#############################################################
# Internal functions.
#############################################################

#' @importFrom Matrix qr qr.coef colSums
#' @importFrom scater nexprs
.per_cluster_normalize <- function(x, curdex, sizes, subset.row, min.mean=1, positive=FALSE, scaling=NULL) 
# Computes the normalization factors _within_ each cluster,
# along with the reference pseudo-cell used for normalization. 
# Written as a separate function so that bplapply operates in the scran namespace.
{
    cur.cells <- length(curdex)
    sizes <- sizes[sizes <= cur.cells]

    vals <- .Call(cxx_subset_and_divide, x, subset.row-1L, curdex-1L, scaling)
    scaling <- vals[[1]]
    exprs <- vals[[2]]
    ave.cell <- vals[[3]] # equivalent to calcAverage().

    # Filtering by mean:
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
    out <- .Call(cxx_forge_system, cur.exprs, ave.cell, sphere - 1L, pool.sizes)
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
    return(rescaling)
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
setGeneric("computeSumFactors", function(x, ...) standardGeneric("computeSumFactors"))

#' @export
setMethod("computeSumFactors", "ANY", .computeSumFactors)

#' @importFrom SummarizedExperiment assay 
#' @importFrom BiocGenerics "sizeFactors<-"
#' @export
setMethod("computeSumFactors", "SingleCellExperiment", function(x, ..., subset.row=NULL, assay.type="counts", get.spikes=FALSE, sf.out=FALSE) 
{ 
    subset.row <- .SCE_subset_genes(subset.row=subset.row, x=x, get.spikes=get.spikes)
    sf <- .computeSumFactors(assay(x, i=assay.type), subset.row=subset.row, ...) 

    if (sf.out) { 
        return(sf) 
    }
    sizeFactors(x) <- sf
    x
})
    
