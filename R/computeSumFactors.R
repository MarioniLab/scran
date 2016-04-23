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

setMethod("computeSumFactors", "ANY", function(x, sizes=c(20, 40, 60, 80, 100), clusters=NULL, ref.clust=NULL, positive=FALSE, errors=FALSE) 
# This contains the function that performs normalization on the summed counts.
# It also provides support for normalization within clusters, and then between
# clusters to make things comparable. It can also switch to linear inverse models
# to ensure that the estimates are non-negative.
#
# written by Aaron Lun
# created 23 November 2015
# last modified 1 December 2015
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

    # Computing the necessary statistics.
    lib.sizes <- colSums(x)
    if (any(lib.sizes < 1e-8)) {
        stop("cells should have non-zero library sizes")
    }
    exprs <- t(t(x)/lib.sizes)
    clust.nf <- clust.profile <- clust.libsize <- list()

    # Computing normalization factors within each cluster first.
    warned.size <- FALSE
    warned.neg <- FALSE

    for (clust in seq_along(indices)) { 
        curdex <- indices[[clust]]
        cur.exprs <- exprs[,curdex,drop=FALSE]
        cur.libs <- lib.sizes[curdex]
        cur.cells <- length(curdex)

        # Checking cluster sizes
        if (any(sizes > cur.cells)) { 
            stop("not enough cells in each cluster for specified 'sizes'") 
        } else if (any(sizes*2L > cur.cells)) {
            if (!warned.size) { warning("number of cells in each cluster should be at least twice that of the largest 'sizes'") }
            warned.size <- TRUE
        }

        # Getting rid of zeros.
        ave.cell <- rowMeans(cur.exprs)
        keep <- ave.cell > .Machine$double.xmin
        use.ave.cell <- ave.cell[keep]
        cur.exprs <- cur.exprs[keep,,drop=FALSE]
        ngenes <- sum(keep)

        # Using our summation approach.
        sphere <- .generateSphere(cur.libs) 
        new.sys <- .create_linear_system(ngenes, cur.cells, cur.exprs, sphere, sizes, use.ave.cell) 
        design <- new.sys$design
        output <- new.sys$output

        # Weighted least-squares (inverse model for positivity).
        if (positive) { 
            design <- as.matrix(design)
            fitted <- limSolve::lsei(A=design, B=output, G=diag(cur.cells), H=numeric(cur.cells), type=2)
            final.nf <- fitted$X
        } else if (errors) {
            design <- as.matrix(design)
            fit <- limma::lmFit(output, design)
            final.nf <- fit$coefficients
            se.est <- fit$stdev.unscaled[1] * fit$sigma # All balanced, so they're all the same.
        } else {
            final.nf <- qr.coef(qr(design), output)
            if (any(final.nf < 0)) { 
                if (!warned.neg) { warning("negative factor estimates, re-run with 'positive=TRUE'") }
                warned.neg <- TRUE
            }
        }

        # Adding per-cluster information.
        clust.nf[[clust]] <- final.nf
        clust.profile[[clust]] <- ave.cell
        clust.libsize[[clust]] <- mean(cur.libs)
    }

    # Adjusting size factors between clusters (using the cluster with the
    # median per-cell library size as the reference, if not specified).
    if (is.null(ref.clust)) {
        clust.libsize <- unlist(clust.libsize)
        ref.col <- which(rank(clust.libsize, ties.method="first")==as.integer(length(clust.libsize)/2)+1L)
    } else {
        ref.col <- which(names(indices)==ref.clust)
        if (length(ref.col)==0L) { 
            stop("'ref.clust' value not in 'clusters'")
        }
    }
    for (clust in seq_along(indices)) { 
        clust.nf[[clust]] <- clust.nf[[clust]] * median(clust.profile[[clust]]/clust.profile[[ref.col]], na.rm=TRUE)
    }
    clust.nf <- unlist(clust.nf)

    # Returning centered size factors, rather than normalization factors.
    final.sf <- rep(NA_integer_, ncells)
    final.sf[unlist(indices)] <- clust.nf
    final.sf <- final.sf * lib.sizes
    
    is.pos <- final.sf > 0 & !is.na(final.sf)
    gm <- exp(mean(log(final.sf[is.pos])))
    final.sf <- final.sf/gm

    if (errors) {
        attr(final.sf, "standard.error") <- se.est
    }
    return(final.sf)
})

setMethod("computeSumFactors", "SCESet", function(x, ..., assay="counts", get.spikes=FALSE) { 
    sf <- computeSumFactors(.getUsedMatrix(x, assay, get.spikes), ...) 
    sizeFactors(x) <- sf
    x
})

.create_linear_system <- function(ngenes, cur.cells, cur.exprs, sphere, sizes, use.ave.cell) {
    sphere <- sphere - 1L # zero-indexing in C++.
    row.dex <- col.dex <- output <- list()
    last.row <- 0L

    for (si in seq_along(sizes)) { 
        out <- .Call(cxx_forge_system, ngenes, cur.cells, cur.exprs, sphere, sizes[si], use.ave.cell)
        if (is.character(out)) { stop(out) }
        row.dex[[si]] <- out[[1]] + last.row
        col.dex[[si]] <- out[[2]]
        output[[si]]<- out[[3]]
        last.row <- last.row + cur.cells
    }
    
    # Adding extra equations to guarantee solvability (downweighted).
    out <- .Call(cxx_forge_system, ngenes, cur.cells, cur.exprs, sphere, 1L, use.ave.cell)
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

