.computeSumFactors <- function(x, sizes=seq(20, 100, 5), clusters=NULL, ref.clust=NULL, 
                               positive=FALSE, errors=FALSE, mean.warn=TRUE, subset.row=NULL) 
# This contains the function that performs normalization on the summed counts.
# It also provides support for normalization within clusters, and then between
# clusters to make things comparable. It can also switch to linear inverse models
# to ensure that the estimates are non-negative.
#
# written by Aaron Lun
# created 23 November 2015
# last modified 23 June 2017 
{
    ncells <- ncol(x)
    if (!is.null(clusters)) {
        if (ncells!=length(clusters)) { 
            stop("'x' ncols is not equal to 'clusters' length")
        }
        is.okay <- !is.na(clusters)
        indices <- split(which(is.okay), clusters[is.okay])
    } else {
        indices <- list(seq_len(ncells))
    }

    # Checking sizes.
    sizes <- sort(as.integer(sizes))
    if (anyDuplicated(sizes)) { 
        stop("'sizes' are not unique") 
    }

    # Checking the subsetting.
    subset.row <- .subset_to_index(subset.row, x, byrow=TRUE)

    # Setting some other values.
    nclusters <- length(indices)
    clust.nf <- clust.profile <- clust.libsizes <- clust.meanlib <- clust.se <- vector("list", nclusters)
    warned.neg <- FALSE
    all.lib <- all.profile <- 0

    # Computing normalization factors within each cluster first.
    for (clust in seq_len(nclusters)) { 
        curdex <- indices[[clust]]
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
        use.ave.cell <- cur.out[[3]]
        ave.cell <- cur.out[[4]]

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
                if (!warned.neg) { warning("encountered negative size factor estimates") }
                warned.neg <- TRUE
            }
            if (errors) {
                warning("errors=TRUE is no longer supported")
            }
        }

        # Adding per-cluster information.
        clust.nf[[clust]] <- final.nf
        clust.profile[[clust]] <- ave.cell
        clust.libsizes[[clust]] <- cur.libs
        clust.meanlib[[clust]] <- mean(cur.libs)
        
        # Storing things to get the average count for warning.
        all.lib <- all.lib + sum(cur.libs)
        all.profile <- all.profile + ave.cell * length(cur.libs)
    }

    # Printing a warning if we see low-abundance genes.
    # This mimics the calcAverage calculation, after adjusting for the mean library size.
    if (mean.warn) {
        mean.profile <- all.profile/ncells 
        mean.lib <- all.lib/ncells
        if (any(mean.profile*mean.lib < 0.1)) {
            warning("low-abundance genes (library size-adjusted average counts below 0.1) detected")
        }
    }

    # Adjusting size factors between clusters (using the cluster with the
    # median per-cell library size as the reference, if not specified).
    if (is.null(ref.clust)) {
        clust.meanlib <- unlist(clust.meanlib)
        ref.col <- which(rank(clust.meanlib, ties.method="first")==as.integer(length(clust.meanlib)/2)+1L)
    } else {
        ref.col <- which(names(indices)==ref.clust)
        if (length(ref.col)==0L) { 
            stop("'ref.clust' value not in 'clusters'")
        }
    }
    clust.nf.scaled <- vector("list", nclusters)
    for (clust in seq_len(nclusters)) { 
        clust.nf.scaled[[clust]] <- clust.nf[[clust]] * median(clust.profile[[clust]]/clust.profile[[ref.col]], na.rm=TRUE)
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

.generateSphere <- function(lib.sizes) 
# This function sorts cells by their library sizes, and generates an ordering vector.
{
    nlibs <- length(lib.sizes)
    o <- order(lib.sizes)
    even <- seq(2,nlibs,2)
    odd <- seq(1,nlibs,2)
    out <- c(o[odd], rev(o[even]))
    c(out, out)
}

LOWWEIGHT <- 0.000001

.create_linear_system <- function(cur.exprs, ave.cell, sphere, pool.sizes) {
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

setGeneric("computeSumFactors", function(x, ...) standardGeneric("computeSumFactors"))

setMethod("computeSumFactors", "ANY", .computeSumFactors)

setMethod("computeSumFactors", "SCESet", function(x, subset.row=NULL, ..., assay="counts", get.spikes=FALSE, sf.out=FALSE) { 
    mat <- assayDataElement(x, assay)
    despiked <- .spike_subset(x, get.spikes)
    
    if (is.null(subset.row)) { 
        subset.row <- despiked 
    } else {
        subset.row <- .subset_to_index(subset.row, mat, byrow=TRUE)
        if (!is.null(despiked)) { 
            subset.row <- intersect(subset.row, which(despiked))
        }
    }

    sf <- .computeSumFactors(mat, subset.row=subset.row, ...) 
    if (sf.out) { 
        return(sf) 
    }
    sizeFactors(x) <- sf
    x
})
    
