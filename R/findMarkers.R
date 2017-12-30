.findMarkers <- function(x, clusters, design=NULL, pval.type=c("any", "all"), direction=c("any", "up", "down"), 
                         min.mean=0.1, subset.row=NULL)
# Uses limma to find the markers that are differentially expressed between clusters,
# given a log-expression matrix and some blocking factors in 'design'.
#
# written by Aaron Lun
# created 22 March 2017
{
    # Creating a design matrix.
    clusters <- as.factor(clusters)
    full.design <- model.matrix(~0 + clusters)
    colnames(full.design) <- clust.vals <- levels(clusters)

    if (!is.null(design)) {
        # Removing terms to avoid linear dependencies on the intercept.
        out <- qr.solve(design, cbind(rep(1, nrow(design))))
        to.drop <- abs(out) > 1e-8
        if (any(to.drop)) {
            design <- design[,-which(to.drop)[1],drop=FALSE]
        }
        full.design <- cbind(full.design, design) # Other linear dependencies will trigger errors in .ranksafe_QR. 
    }

    # Other variables to be initialized.
    pval.type <- match.arg(pval.type) 
    direction <- match.arg(direction)  
    subset.row <- .subset_to_index(subset.row, x, byrow=TRUE)

    # Estimating the parameters.
    fit <- .fit_lm_internal(x=x, subset.row=subset.row, clusters=clusters, design=full.design)
    coefficients <- fit$coefficients
    means <- fit$means
    sigma2 <- fit$sigma2

    # Performing the EB shrinkage in two chunks if min.mean is requested.
    # This ensures that discreteness at low means is quarantined from the high abundances.
    higher <- means >= min.mean
    if (is.null(min.mean) || all(higher) || !any(higher)) { 
        eb.out <- .perform_eb_shrinkage(sigma2, covariate=means, design=full.design)
    } else {
        high.out <- .perform_eb_shrinkage(sigma2[higher], covariate=means[higher], design=full.design)
        low.out <- .perform_eb_shrinkage(sigma2[!higher], covariate=means[!higher], design=full.design)
        eb.out <- mapply(low.out, high.out, FUN=function(low, high) {
            val <- numeric(length(means))
            val[higher] <- high
            val[!higher] <- low
            val            
        }, SIMPLIFY=FALSE)
    }

    # Doing a dummy fit, to avoid having to manually calculate standard errors.
    lfit <- lmFit(rbind(seq_len(nrow(full.design))), full.design)

    output <- vector("list", length(clust.vals))
    names(output) <- clust.vals  
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
            cur.t <- cur.lfc/(lfit2$stdev.unscaled[con]*sqrt(eb.out$var.post))
            cur.p <- 2 * pt(-abs(cur.t), df=eb.out$df.total)

            if (direction=="up") { 
                cur.p <- ifelse(cur.lfc > 0, cur.p/2, 1-cur.p/2)
            } else if (direction=="down") {
                cur.p <- ifelse(cur.lfc < 0, cur.p/2, 1-cur.p/2)
            }
            all.p[,con] <- cur.p
        }
        colnames(all.lfc) <- paste0("logFC.", targets)

        # Collating p-values.
        pval <- .combine_pvalues(all.p, pval.type)
        rank.out <- .rank_top_genes(all.p)
        min.rank <- rank.out$rank
        min.p <- rank.out$value

        # Producing the output object.
        marker.set <- data.frame(Top=min.rank, Gene=rownames(x)[subset.row], 
                                 FDR=p.adjust(pval, method="BH"), all.lfc,
                                 stringsAsFactors=FALSE, check.names=FALSE)
        marker.set <- marker.set[order(marker.set$Top, min.p),]
        rownames(marker.set) <- NULL
        output[[host]] <- marker.set
    }

    return(output)
}

###########################################################
# Internal functions.
###########################################################

.fit_lm_internal <- function(x, subset.row, clusters, design) 
# Fit a linear model and get coefficients (need to account for the fact that the matrix is pivoted).
# Alternatively, do this blockwise to avoid having to do the QtY multiplication and back-solving.
{
    if (nlevels(clusters) != ncol(design)) {
        QR <- .ranksafe_qr(design)
        if (nrow(design) <= ncol(design)) {
            stop("no residual d.f. in design matrix for variance estimation") 
        }
        stats <- .Call(cxx_fit_linear_model, QR$qr, QR$qraux, x, subset.row - 1L, TRUE)
        coefficients <- stats[[1]][order(QR$pivot),,drop=FALSE]
        means <- stats[[2]]
        sigma2 <- stats[[3]]
    } else {
        # Assumes ordering of levels in 'clusters' is the same as that in 'design'.
        by.block <- split(seq_len(ncol(x))-1L, clusters, drop=TRUE)
        group.size <- lengths(by.block)
        resid.df <- group.size - 1L
        if (all(resid.df<=0L)){ 
            stop("no residual d.f. in any level of 'clusters' for variance estimation")
        }

        # Calculating the statistics for each block, and adding them across blocks
        # (avoid NA variances when some blocks have no residual d.f.).
        stats <- .Call(cxx_fit_oneway, by.block, x, subset.row - 1L)
        coefficients <- stats[[1]]
        means <- drop(coefficients %*% group.size)/sum(group.size)
        keep <- resid.df > 0
        sigma2 <- drop(stats[[2]][,keep,drop=FALSE] %*% resid.df[keep])/sum(resid.df[keep])
        coefficients <- t(coefficients) # for consistency.
    }
    return(list(coefficients=coefficients, means=means, sigma2=sigma2))
}

.perform_eb_shrinkage <- function(sigma2, covariate, design)
# EB shrinkage, limma-style.
{
    df.residual <- rep(nrow(design) - ncol(design), length(covariate))
    eb.out <- squeezeVar(sigma2, df=df.residual, robust=TRUE, covariate=covariate)
    df.total <- df.residual + eb.out$df.prior
    df.pooled <- sum(df.residual, na.rm = TRUE)
    df.total <- pmin(df.total, df.pooled)
    return(list(var.post=eb.out$var.post, df.total=df.total))
}

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
    ngenes <- nrow(all.p)
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


