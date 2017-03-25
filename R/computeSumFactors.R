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

setGeneric("computeSumFactors", function(x, ...) standardGeneric("computeSumFactors"))

setMethod("computeSumFactors", "matrix", function(x, sizes=c(20, 40, 60, 80, 100), clusters=NULL, ref.clust=NULL, positive=FALSE, errors=FALSE, subset.row=NULL) 
# This contains the function that performs normalization on the summed counts.
# It also provides support for normalization within clusters, and then between
# clusters to make things comparable. It can also switch to linear inverse models
# to ensure that the estimates are non-negative.
#
# written by Aaron Lun
# created 23 November 2015
# last modified 22 March 2017
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
    sizes <- as.integer(sizes)
    if (anyDuplicated(sizes)) { 
        stop("'sizes' is not unique") 
    }

    # Checking the subsetting.
    subset.row <- .subset_to_index(subset.row, x, byrow=TRUE)

    # Setting some other values.
    nclusters <- length(indices)
    clust.nf <- clust.profile <- clust.libsizes <- clust.meanlib <- vector("list", nclusters)
    warned.neg <- FALSE

    # Computing normalization factors within each cluster first.
    for (clust in seq_len(nclusters)) { 
        curdex <- indices[[clust]]
        cur.out <- .Call(cxx_subset_and_divide, x, subset.row-1L, curdex-1L) 
        if (is.character(cur.out)) { stop(cur.out) }
        cur.libs <- cur.out[[1]]
        cur.exprs <- cur.out[[2]]       
        cur.cells <- length(curdex)

        # Checking cluster sizes
        if (any(sizes > cur.cells)) { 
            stop("not enough cells in each cluster for specified 'sizes'") 
        } 

        # Getting rid of zeros.
        ave.cell <- rowMeans(cur.exprs)
        keep <- ave.cell > .Machine$double.xmin
        use.ave.cell <- ave.cell[keep]
        cur.exprs <- cur.exprs[keep,,drop=FALSE]

        # Using our summation approach.
        sphere <- .generateSphere(cur.libs) 
        new.sys <- .create_linear_system(cur.exprs, sphere, sizes, use.ave.cell) 
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
                if (!warned.neg) { warning("negative factor estimates, re-run with 'positive=TRUE'") }
                warned.neg <- TRUE
            }

            if (errors) {
                # Our "observations" here _are_ our size factors, so variance refers to that of the size factors.
                # Don't compute the standard error of the coefficients, as that isn't the relevant value here.
                sigma2 <- mean(qr.qty(QR, output)[-seq_len(ncol(design))]^2)
                se.est <- sqrt(sigma2)
            }
        }

        # Adding per-cluster information.
        clust.nf[[clust]] <- final.nf
        clust.profile[[clust]] <- ave.cell
        clust.libsizes[[clust]] <- cur.libs
        clust.meanlib[[clust]] <- mean(cur.libs)
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

    if (errors) {
        # Calculating all the changes to get to this point, and scaling the standard error by them.
        clust.nf <- unlist(clust.nf)
        clust.nf[indices] <- clust.nf
        attr(final.sf, "standard.error") <- se.est * final.sf/clust.nf
    }
    return(final.sf)
})

setMethod("computeSumFactors", "SCESet", function(x, subset.row=NULL, ..., assay="counts", get.spikes=FALSE, sf.out=FALSE) { 
    if (is.null(subset.row)) { 
        subset.row <- .spikeSubset(x, get.spikes)
    }
    sf <- computeSumFactors(assayDataElement(x, assay), subset.row=subset.row, ...) 
    if (sf.out) { 
        return(sf) 
    }
    sizeFactors(x) <- sf
    x
})

.create_linear_system <- function(cur.exprs, sphere, sizes, use.ave.cell) {
    sphere <- sphere - 1L # zero-indexing in C++.
    row.dex <- col.dex <- output <- vector("list", length(sizes))
    last.row <- 0L
    cur.cells <- ncol(cur.exprs)

    for (si in seq_along(sizes)) { 
        out <- .Call(cxx_forge_system, cur.exprs, sphere, sizes[si], use.ave.cell)
        if (is.character(out)) { stop(out) }
        row.dex[[si]] <- out[[1]] + last.row
        col.dex[[si]] <- out[[2]]
        output[[si]]<- out[[3]]
        last.row <- last.row + cur.cells
    }
    
    # Adding extra equations to guarantee solvability (downweighted).
    out <- .Call(cxx_forge_system, cur.exprs, sphere, 1L, use.ave.cell)
    if (is.character(out)) { stop(out) }
    si <- length(row.dex) + 1L
    row.dex[[si]] <- out[[1]] + last.row
    col.dex[[si]] <- out[[2]]
    output[[si]] <- out[[3]]

    # Weighting the system.
    LOWWEIGHT <- 0.000001
    output[[si]] <- output[[si]] * sqrt(LOWWEIGHT)
    eqn.values <- rep(rep(c(1, sqrt(LOWWEIGHT)), c(si-1L, 1L)), lengths(row.dex))

    # Constructing a sparse matrix.
    row.dex <- unlist(row.dex)
    col.dex <- unlist(col.dex)
    eqn.values <- unlist(eqn.values)
    output <- unlist(output)
    design <- sparseMatrix(i=row.dex + 1L, j=col.dex + 1L, x=eqn.values, dims=c(length(output), cur.cells))

    return(list(design=design, output=output))
}

