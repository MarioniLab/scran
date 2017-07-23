.findMarkers <- function(x, clusters, design=NULL, pval.type=c("any", "all"), direction=c("any", "up", "down"), 
                         min.mean=0.1, subset.row=NULL)
# Uses limma to find the markers that are differentially expressed between clusters,
# given a log-expression matrix and some blocking factors in 'design'.
#
# written by Aaron Lun
# created 22 March 2017
# last modified 22 July 2017    
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
        full.design <- cbind(full.design, design) # Other linear dependencies will trigger warnings.
    }
 
    pval.type <- match.arg(pval.type) 
    direction <- match.arg(direction)  
  
    # Fit a linear model with a one-way layout (need to account for the fact that the matrix is pivoted).
    subset.row <- .subset_to_index(subset.row, x, byrow=TRUE)
    QR <- .ranksafe_qr(full.design)
    stats <- .Call(cxx_fit_linear_model, QR$qr, QR$qraux, x, subset.row - 1L, TRUE)
    coefficients <- stats[[1]][order(QR$pivot),,drop=FALSE]
    means <- stats[[2]]
    sigma2 <- stats[[3]]

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
            largest <- (max.col(all.p) - 1) * ngenes + seq_len(ngenes)
            pval <- all.p[largest]
        }

        # Collating minimum ranks.
        min.rank <- rep(ngenes, ngenes)
        min.p <- rep(1, ngenes)
        for (con in seq_len(ncon)) { 
            cur.p <- all.p[,con]
            cur.rank <- rank(cur.p, ties.method="first")
            min.rank <- pmin(min.rank, cur.rank)
            min.p <- pmin(min.p, cur.p)
        }

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

setGeneric("findMarkers", function(x, ...) standardGeneric("findMarkers"))

setMethod("findMarkers", "matrix", .findMarkers)

setMethod("findMarkers", "SCESet", function(x, ..., subset.row=NULL, assay="exprs", get.spikes=FALSE) {
    if (is.null(subset.row)) { subset.row <- .spike_subset(x, get.spikes) }
    .findMarkers(assayDataElement(x, assay), ..., subset.row=subset.row)
})                                 


