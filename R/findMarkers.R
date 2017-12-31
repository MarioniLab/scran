.findMarkers <- function(x, clusters, block=NULL, design=NULL, pval.type=c("any", "all"), direction=c("any", "up", "down"), subset.row=NULL)
# Uses limma to find the markers that are differentially expressed between clusters,
# given a log-expression matrix and some blocking factors in 'design'.
#
# written by Aaron Lun
# created 22 March 2017
{
    ncells <- ncol(x)
    clusters <- as.factor(clusters)
    if (length(clusters)!=ncells) {
        stop("length of 'clusters' does not equal 'ncol(x)'")
    }
    pval.type <- match.arg(pval.type) 
    direction <- match.arg(direction)  
    subset.row <- .subset_to_index(subset.row, x, byrow=TRUE)

    # Estimating the parameters.
    if (!is.null(block) || is.null(design)) {
        fit <- .test_block_internal(x, subset.row, clusters, block, direction)
    } else {
        fit <- .fit_lm_internal(x, subset.row, clusters, design, direction) 
    }

    # Figuring out how to rank within each DataFrame.
    all.p <- fit$p.value
    all.lfc <- fit$logFC
    output <- vector("list", length(all.p))
    names(output) <- names(all.p)

    for (host in names(all.p)) {
        pval <- .combine_pvalues(all.p[[host]], pval.type)
        lfc <- all.lfc[[host]]
        colnames(lfc) <- sprintf("logFC.%s", colnames(lfc))

        if (pval.type=="any") {
            rank.out <- .rank_top_genes(all.p[[host]])
            min.rank <- rank.out$rank
            min.p <- rank.out$value
            o <- order(min.rank, min.p)
            preamble <- DataFrame(Top=min.rank)
        } else {
            o <- order(pval)
            preamble <- DataFrame(Worst=pval)
        }

        # Producing the output object.
        gene.names <- rownames(x)[subset.row]
        if (is.null(gene.names)) { 
            gene.names <- subset.row
        }
        marker.set <- DataFrame(preamble, FDR=p.adjust(pval, method="BH"), lfc,
                                check.names=FALSE, row.names=gene.names)
        marker.set <- marker.set[o,,drop=FALSE]
        output[[host]] <- marker.set
    }

    return(output)
}

###########################################################
# Internal functions (blocking)
###########################################################

.test_block_internal <- function(x, subset.row, clusters, block, direction) 
# This looks at every level of the blocking factor and performs
# t-tests between pairs of clusters within each blocking level.
{
    ncells <- ncol(x)
    if (is.null(block)) {
        by.block <- list(`1`=seq_len(ncells))
    } else {
        if (length(block)!=ncells) {
            stop("length of 'block' does not equal 'ncol(x)'")
        } 
        by.block <- split(seq_len(ncells), block)
    }
        
    # Setting up block-level outputs. 
    out.means <- vector("list", length(by.block))
    out.s2 <- out.n <- out.means
    clust.vals <- levels(clusters)

    for (b in seq_along(by.block)) {
        curblock <- by.block[[b]]       
        by.cluster <- split(curblock - 1L, clusters[curblock], drop=FALSE)
        stopifnot(identical(names(by.cluster), clust.vals))
     
        # Calculating the statistics for each block.
        stats <- .Call(cxx_fit_oneway, by.cluster, x, subset.row - 1L)
        means <- stats[[1]]
        sigma2 <- stats[[2]]
        colnames(means) <- colnames(sigma2) <- clust.vals

        out.means[[b]] <- means
        out.s2[[b]] <- sigma2
        out.n[[b]] <- lengths(by.cluster)
    }

    # Setting up output containers.
    out.lfc <- vector("list", length(clust.vals))
    names(out.lfc) <- clust.vals  
    out.p <- out.lfc
    ngenes <- length(subset.row)
    for (i in seq_along(clust.vals)) {
        targets <- clust.vals[-i]
        collected <- matrix(0, ngenes, length(targets))
        colnames(collected) <- targets
        host <- clust.vals[i]
        out.p[[host]] <- out.lfc[[host]] <- collected
    }

    # Running through all pairs of comparisons.
    nblocks <- length(by.block)
    for (i in seq_along(clust.vals)) {
        host <- clust.vals[i]
        targets <- clust.vals[seq_len(i-1L)]

        for (target in targets) {
            all.lfc <- matrix(0, ngenes, nblocks)
            all.left <- all.right <- all.weight <- all.lfc
            all.rss <- all.df <- numeric(ngenes)

            # Performing the same pairwise t-test across all blocks.
            for (b in seq_len(nblocks)) { 
                cur.lfc <- out.means[[b]][,host] - out.means[[b]][,target]
                t.out <- .get_t_test_values(host.s2=out.s2[[b]][,host], target.s2=out.s2[[b]][,target],
                                            host.n=out.n[[b]][host], target.n=out.n[[b]][target])
                cur.err <- t.out$err
                cur.df <- t.out$test.df
                cur.t <- cur.lfc/sqrt(cur.err)

                all.lfc[,b] <- cur.lfc
                all.weight[,b] <- 1/cur.err
                all.left[,b] <- pt(cur.t, df=cur.df, lower.tail=TRUE, log.p=TRUE)
                all.right[,b] <- pt(cur.t, df=cur.df, lower.tail=FALSE, log.p=TRUE)

                all.rss <- all.rss + t.out$rss 
                all.df <- all.df + t.out$residual.df
            }

            # Combining the p-values and log-fold changes across blocks.
            com.p <- .combine_test_statistics(all.left, all.right, all.weight, direction)
            com.lfc <- .combine_lfc(all.lfc, all.weight, all.rss/all.df)
            out.p[[host]][,target] <- com.p
            out.lfc[[host]][,target] <- com.lfc

            # Filling in the other values due to symmetry. 
            out.lfc[[target]][,host] <- -com.lfc
            if (direction=="any") { 
                out.p[[target]][,host] <- com.p
            } else {
                out.p[[target]][,host] <- .combine_test_statistics(all.right, all.left, all.weight, direction)
            }
        }
    }
    return(list(p.value=out.p, logFC=out.lfc))    
}

.get_t_test_values <- function(host.s2, target.s2, host.n, target.n) 
# Computes the error variance of the log-fold change and the 
# degrees of freedom for the t-distribution.
{
    host.df <- max(0L, host.n - 1L)
    target.df <- max(0L, target.n - 1L)

    if (host.df > 0L && target.df > 0L) { 
        # Perform Welch's t-test here.
        host.err <- host.s2/host.n
        target.err <- target.s2/target.n
        cur.err <- host.err + target.err
        cur.df <- cur.err^2 / (host.err^2/host.df + target.err^2/target.df)
       
    } else {
        if (host.n==0L || target.n==0L) { 
            cur.err <- rep(NA_real_, length(host.s2))
            cur.df <- NA_real_

        } else {
            # Try to perform Student's t-test, if possible.
            if (host.df==0L) {
                cur.err <- target.s2 * (1 + 1/target.n)
                cur.df <- target.df
            } else {
                cur.err <- host.s2 * (1 + 1/host.n)
                cur.df <- host.df
            }
        }
    }
  
    # Avoid unlikely but potential problems with discreteness. 
    cur.err <- pmax(cur.err, 1e-8) 

    # Also computing some general variance statistics.
    residual.df <- host.df + target.df
    rss <- ifelse(is.na(host.s2), 0, host.s2) * host.df + ifelse(is.na(target.s2), 0, target.s2) * target.df 
    return(list(err=cur.err, test.df=cur.df, rss=rss, residual.df=residual.df))
}

.combine_lfc <- function(lfc, weights, s2) 
# This computes the weighted average log-fold change across all batches.
# Weighting is usually based on the squared error of the log-fold change,
# so some finesse is required to use information from batches with n=1
# for all groups.
{
    if (length(weights) && is.na(max(weights))) {
        s2[is.na(s2)] <- 1 # effectively ensures equal weights, if no one has any residual d.f.
        for (i in seq_len(ncol(weights))) {
            curw <- weights[,i]
            lost <- is.na(curw) & !is.na(lfc[,i])
            curw[lost] <- 1/(s2[lost] * 2) # SE^2 when n1=n2=1 with known SD (set to average across blocks).
            weights[,i] <- curw
        }
    }
    rowSums(weights*lfc, na.rm=TRUE)/rowSums(weights, na.rm=TRUE)
}

.combine_test_statistics <- function(left, right, weights, direction) 
# This uses the inverse errors as the weights to combine the 
# various statistics. It uses Stouffer's method to combine 
# the p-values on either side.
{
    sum.weight <- rowSums(weights, na.rm=TRUE)
    if (direction!="up") { 
        Z.left <- qnorm(left, log.p=TRUE)
        final.left <- rowSums(Z.left*weights, na.rm=TRUE)/sum.weight
        p.out <- p.left <- pnorm(final.left)
    }

    if (direction!="down") { 
        Z.right <- qnorm(right, log.p=TRUE)
        final.right <- rowSums(Z.right*weights, na.rm=TRUE)/sum.weight
        p.out <- p.right <- pnorm(final.right)
    }
    
    if (direction=="any") {
        p.out <- pmin(p.left, p.right)*2
    }
    return(p.out)
}

###########################################################
# Internal functions (linear modelling)
###########################################################

.fit_lm_internal <- function(x, subset.row, clusters, design, direction) 
# This fits a linear model to each gene and performs a t-test for 
# differential expression between clusters.
{
    full.design <- model.matrix(~0 + clusters)
    clust.vals <- levels(clusters)
    colnames(full.design) <- clust.vals 

    # Removing terms to avoid linear dependencies on the intercept.
    if (nrow(design)!=ncol(x)) { 
        stop("'nrow(design)' is not equal to 'ncol(x)'")
    }
    out <- qr.solve(design, cbind(rep(1, nrow(design))))
    to.drop <- abs(out) > 1e-8
    if (any(to.drop)) {
        design <- design[,-which(to.drop)[1],drop=FALSE]
    }
    full.design <- cbind(full.design, design) # Other linear dependencies will trigger errors in .ranksafe_QR. 

    # Getting coefficient estimates (need to undo column pivoting to get the actual estimates).
    QR <- .ranksafe_qr(full.design)
    resid.df <- nrow(full.design) - ncol(full.design)
    if (resid.df <= 0L) {
        stop("no residual d.f. in design matrix for variance estimation") 
    }
    stats <- .Call(cxx_fit_linear_model, QR$qr, QR$qraux, x, subset.row - 1L, TRUE)
    coefficients <- stats[[1]]
    coefficients[QR$pivot,] <- coefficients
    # means <- stats[[2]] # Not actually needed here anymore.
    sigma2 <- stats[[3]]

    # Setting up output containers.
    out.lfc <- vector("list", length(clust.vals))
    names(out.lfc) <- clust.vals  
    out.p <- out.lfc

    # Doing a dummy fit, to avoid having to manually calculate standard errors.
    lfit <- lmFit(rbind(seq_len(nrow(full.design))), full.design)
    for (h in seq_along(clust.vals)) { 
        host <- clust.vals[h]
        not.h <- seq_along(clust.vals)[-h]
        targets <- clust.vals[not.h]
        
        # Computing standard errors (via limma).
        con <- matrix(0, ncol(full.design), length(clust.vals))
        diag(con) <- -1
        con[h,] <- 1
        con <- con[,not.h,drop=FALSE]
        lfit2 <- contrasts.fit(lfit, con)

        # Computing log-fold changes, t-statistics and p-values on a per-contrast basis.
        # This _could_ be vectorised, but it's less confusing to do it like this,
        # and there's not much speed gain to be had from vectorizing over contrasts.
        ngenes <- length(subset.row)
        ncon <- length(not.h)
        all.lfc <- all.p <- matrix(0, ngenes, ncon)
        ref.coef <- coefficients[h,]

        for (con in seq_len(ncon)) { 
            cur.lfc <- ref.coef - coefficients[not.h[con],]
            all.lfc[,con] <- cur.lfc
            cur.t <- cur.lfc/(lfit2$stdev.unscaled[con]*sqrt(sigma2))
            cur.p <- 2 * pt(-abs(cur.t), df=resid.df)

            if (direction=="up") { 
                cur.p <- ifelse(cur.lfc > 0, cur.p/2, 1-cur.p/2)
            } else if (direction=="down") {
                cur.p <- ifelse(cur.lfc < 0, cur.p/2, 1-cur.p/2)
            }
            all.p[,con] <- cur.p
        }

        colnames(all.p) <- colnames(all.lfc) <- targets
        out.p[[host]] <- all.p
        out.lfc[[host]] <- all.lfc
    }
    return(list(p.value=out.p, logFC=out.lfc))
}

###########################################################
# Internal functions (other)
###########################################################

.rank_top_genes <- function(metrics) 
# This computes the rank and the minimum metric for each gene.
{
    ngenes <- nrow(metrics)
    ncon <- ncol(metrics)
    min.rank <- rep(ngenes, ngenes)
    min.val <- rep(1, ngenes)

    for (con in seq_len(ncon)) { 
        cur.val <- metrics[,con]
        cur.rank <- rank(cur.val, ties.method="first")
        min.rank <- pmin(min.rank, cur.rank)
        min.val <- pmin(min.val, cur.val)
    }
    
    return(list(rank=min.rank, value=min.val))
}

.combine_pvalues <- function(all.p, pval.type) { 
    # Getting rid of NA's.
    ngenes <- nrow(all.p)
    if (ngenes) {
        discard <- is.na(all.p[1,])
        all.p <- all.p[,!discard,drop=FALSE]
    }
    ncon <- ncol(all.p)

    if (pval.type=="any") { 
        # Computing Simes' p-value in a fully vectorised manner.
        gene.id <- rep(seq_len(ngenes), ncon)
        o <- order(gene.id, all.p)
        penalty <- rep(ncon/seq_len(ncon), ngenes) 
        com.p <- matrix(all.p[o]*penalty, ngenes, ncon, byrow=TRUE)
        smallest <- (max.col(-com.p) - 1) * ngenes + seq_len(ngenes)
        pval <- com.p[smallest]
    } else {
        # Computing the IUT p-value.
        pval <- all.p[.find_largest_index(all.p)]
    }

    return(pval)
}

.find_largest_index <- function(metrics) {
    ngenes <- nrow(metrics)
    (max.col(metrics) - 1L) * ngenes + seq_len(ngenes)
}

###########################################################
# S4 method definitions
###########################################################

setGeneric("findMarkers", function(x, ...) standardGeneric("findMarkers"))

setMethod("findMarkers", "ANY", .findMarkers)

setMethod("findMarkers", "SingleCellExperiment", 
          function(x, ..., subset.row=NULL, assay.type="logcounts", get.spikes=FALSE) {

    subset.row <- .SCE_subset_genes(subset.row=subset.row, x=x, get.spikes=get.spikes)
    .findMarkers(assay(x, i=assay.type), ..., subset.row=subset.row)
})                                 


