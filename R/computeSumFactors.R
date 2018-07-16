#' @importFrom BiocParallel bplapply SerialParam
.computeSumFactors <- function(x, sizes=seq(20, 100, 5), clusters=NULL, ref.clust=NULL, max.cluster.size=3000, 
                               positive=FALSE, errors=FALSE, min.mean=1, subset.row=NULL, BPPARAM=SerialParam())
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

    if (length(indices)==0L) {
        # To ensure that the empty cluster error triggers.
        indices <- list(integer(0))
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
        min.mean=min.mean, positive=positive, errors=errors, BPPARAM=BPPARAM)

    clust.nf <- lapply(all.norm, "[[", i="final.nf")
    clust.profile <- lapply(all.norm, "[[", i="ave.cell")
    clust.libsizes <- lapply(all.norm, "[[", i="cur.libs")
    clust.meanlib <- vapply(all.norm, FUN="[[", i="mean.lib", FUN.VALUE=0)

    # Adjusting size factors between clusters.
    rescaling.factors <- .rescale_clusters(clust.profile, clust.meanlib, ref.clust=ref.clust, 
                                           min.mean=min.mean, clust.names=names(indices))
    clust.nf.scaled <- vector("list", nclusters)
    for (clust in seq_len(nclusters)) { 
        clust.nf.scaled[[clust]] <- clust.nf[[clust]] * rescaling.factors[[clust]]
    }
    clust.nf.scaled <- unlist(clust.nf.scaled)

    # Returning centered size factors, rather than normalization factors.
    clust.sf <- clust.nf.scaled * unlist(clust.libsizes) 
    final.sf <- rep(NA_integer_, ncells)
    indices <- unlist(indices)
    final.sf[indices] <- clust.sf
    
    is.pos <- final.sf > 0 & !is.na(final.sf)
    final.sf <- final.sf/mean(final.sf[is.pos])
    return(final.sf)
}

#############################################################
# Internal functions.
#############################################################

#' @importFrom Matrix qr qr.coef
.per_cluster_normalize <- function(x, curdex, sizes, subset.row, min.mean=1, positive=FALSE, errors=FALSE) 
# Computes the normalization factors _within_ each cluster,
# along with the reference pseudo-cell used for normalization. 
# Written as a separate function so that bplapply operates in the scran namespace.
{
    cur.cells <- length(curdex)
    cur.sizes <- sizes
    if (any(cur.sizes > cur.cells)) { 
        cur.sizes <- cur.sizes[cur.sizes <= cur.cells]
        if (length(cur.sizes)) { 
            warning("not enough cells in at least one cluster for some 'sizes'")
        } else {
            stop("not enough cells in at least one cluster for any 'sizes'")
        }
    } 

    cur.out <- .Call(cxx_subset_and_divide, x, subset.row-1L, curdex-1L) 
    cur.libs <- cur.out[[1]]
    exprs <- cur.out[[2]]
    ave.cell <- cur.out[[3]]

    # Filtering by mean (easier to do it here in R, rather than C++).
    mean.lib <- mean(cur.libs)
    high.ave <- min.mean/mean.lib <= ave.cell # mimics calcAverage
    if (!all(high.ave)) { 
        exprs <- exprs[high.ave,,drop=FALSE]
        use.ave.cell <- ave.cell[high.ave]
    } else {
        use.ave.cell <- ave.cell
    }

    # Using our summation approach.
    sphere <- .generateSphere(cur.libs)
    new.sys <- .create_linear_system(exprs, use.ave.cell, sphere, cur.sizes) 
    design <- new.sys$design
    output <- new.sys$output

    # Weighted least-squares (inverse model for positivity).
    if (positive) { 
        design <- as.matrix(design)
        fitted <- limSolve::lsei(A=design, B=output, G=diag(cur.cells), H=numeric(cur.cells), type=2)
        final.nf <- fitted$X
    } else {
        QR <- qr(design)
        final.nf <- qr.coef(QR, output)
        if (any(final.nf < 0)) { 
            warning("encountered negative size factor estimates")
        }
        if (errors) {
            warning("errors=TRUE is no longer supported")
        }
    }

    return(list(final.nf=final.nf, ave.cell=ave.cell, cur.libs=cur.libs, mean.lib=mean.lib))
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
    sphere <- sphere - 1L # zero-indexing in C++.

    nsizes <- length(pool.sizes)
    row.dex <- col.dex <- output <- vector("list", 2L)
    cur.cells <- ncol(cur.exprs)

    # Creating the linear system with the requested pool sizes.
    out <- .Call(cxx_forge_system, cur.exprs, ave.cell, sphere, pool.sizes)
    row.dex[[1]] <- out[[1]] 
    col.dex[[1]] <- out[[2]]
    output[[1]]<- out[[3]]
    
    # Adding extra equations to guarantee solvability (downweighted).
    out <- .Call(cxx_forge_system, cur.exprs, ave.cell, sphere, 1L)
    row.dex[[2]] <- out[[1]] + cur.cells * nsizes
    col.dex[[2]] <- out[[2]]
    output[[2]] <- out[[3]] * sqrt(LOWWEIGHT)

    # Setting up the entries of the LHS matrix.
    eqn.values <- rep(c(1, sqrt(LOWWEIGHT)), lengths(row.dex))

    # Constructing a sparse matrix.
    row.dex <- unlist(row.dex)
    col.dex <- unlist(col.dex)
    output <- unlist(output)
    design <- sparseMatrix(i=row.dex + 1L, j=col.dex + 1L, x=eqn.values, dims=c(length(output), cur.cells))

    return(list(design=design, output=output))
}

#' @importFrom stats median
.rescale_clusters <- function(mean.prof, mean.lib, ref.clust, min.mean, clust.names) 
# Chooses a cluster as a reference and rescales all other clusters to the reference,
# based on the 'normalization factors' computed between pseudo-cells.
{
    # Picking the cluster with the median library size as the reference.
    if (is.null(ref.clust)) {
        ref.col <- which(rank(mean.lib, ties.method="first")==as.integer(length(mean.lib)/2)+1L)
    } else {
        ref.col <- which(clust.names==ref.clust)
        if (length(ref.col)==0L) { 
            stop("'ref.clust' value not in 'clusters'")
        }
    }

    nclusters <- length(mean.prof)
    rescaling <- vector("list", nclusters)
    for (clust in seq_len(nclusters)) { 
        ref.prof <- mean.prof[[ref.col]]
        cur.prof <- mean.prof[[clust]] 

        # Filtering based on the mean of the per-cluster means (requires scaling for the library size).
        # Effectively equivalent to 'calcAverage(cbind(ref.ave.count, cur.ave.count))' where the averages
        # are themselves equivalent to 'calcAverage()' across all cells in each cluster.
        ref.libsize <- mean.lib[[ref.col]]
        cur.libsize <- mean.lib[[clust]]
        to.use <- (cur.prof + ref.prof)/2 * (cur.libsize + ref.libsize)/2 >= min.mean
        if (!all(to.use)) { 
            cur.prof <- cur.prof[to.use]
            ref.prof <- ref.prof[to.use]
        } 

        # Adjusting for systematic differences between clusters.
        rescaling[[clust]] <- median(cur.prof/ref.prof, na.rm=TRUE)
    }
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
setMethod("computeSumFactors", "SingleCellExperiment", 
          function(x, ..., min.mean=1, subset.row=NULL, 
                   assay.type="counts", get.spikes=FALSE, sf.out=FALSE) { 
 
    if (is.null(subset.row) && !is.null(min.mean)) {
        no.warn <- TRUE # avoid triggering the message of the subsetting was not originally specified
    }

    subset.row <- .SCE_subset_genes(subset.row=subset.row, x=x, get.spikes=get.spikes)
    sf <- .computeSumFactors(assay(x, i=assay.type), subset.row=subset.row, 
                             min.mean=min.mean, ...) 

    if (sf.out) { 
        return(sf) 
    }
    sizeFactors(x) <- sf
    x
})
    
