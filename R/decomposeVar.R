.decompose_var <- function(x, fit, block=NA, design=NA, subset.row=NULL, ...)
# Computes the biological variability of the log-CPMs by subtracting the
# inferred technical variance from the total variance.
#
# written by Aaron Lun
# created 21 January 2016 
{   
    if (missing(x)) {
        # Extracting statistics from an existing fit.
        lmeans <- fit$mean
        lvar <- fit$var
        resid.df <- fit$resid.df

        block <- fit$block
        design <- fit$design
        if (!is.null(block)) { 
            ncells <- length(block)
        } else if (!is.null(design)) {
            ncells <- nrow(design)
        } else {
            ncells <- unname(resid.df + 1L)[1]
        }
           
    } else {
        # Computing statistics fresh.
        if (length(block)==1L && is.na(block)) {
            block <- fit$block
        }
        if (length(design)==1L && is.na(design)) {
            design <- fit$design
        }

        ncells <- ncol(x)
        stats.out <- .get_var_stats(x, block=block, design=design, subset.row=subset.row)
        lmeans <- stats.out$means
        lvar <- stats.out$var
        resid.df <- stats.out$resid.df
    }

    if (!is.null(block)) {
        gnames <- rownames(lmeans)

        # Need to do some work to calculate the technical component separately for each block.
        block.tech.var <- fit$trend(as.vector(lmeans))
        block.pval <- testVar(as.vector(lvar), null=block.tech.var, df=as.vector(resid.df), second.df=fit$df2, log.p=TRUE, ...)
        dim(block.pval) <- dim(block.tech.var) <- dim(lmeans)

        # Aggregating variances together, using the weighted mean of residual d.f. 
        # (we set na.rm=TRUE to ignore blocks where no resid. d.f. are available).
        total.resid.df <- rowSums(resid.df)
        lvar <- rowSums(lvar * resid.df, na.rm=TRUE)/total.resid.df
        tech.var <- rowSums(block.tech.var * resid.df)/total.resid.df
        
        # Aggregating means together (effectively rowMeans(x), without having to recompute it from 'x').
        num.per.block <- table(block)
        lmeans <- lmeans %*% num.per.block[colnames(lmeans)]/length(block)

        # Combining p-values using Stouffer's method (independent observations).
        Z <- rowSums(qnorm(block.pval, log.p=TRUE) * resid.df)/total.resid.df
        pval <- pnorm(Z)
        resid.df <- total.resid.df

    } else {
        gnames <- names(lmeans)
        tech.var <- fit$trend(lmeans)
        pval <- testVar(total=lvar, null=tech.var, df=resid.df, second.df=fit$df2, ...)
    }
    
    pval <- as.numeric(pval)
    lvar <- as.numeric(lvar)
    tech.var <- as.numeric(tech.var)
    bio.var <- lvar - tech.var
    out <- DataFrame(mean=as.numeric(lmeans), 
                     total=lvar, bio=bio.var, tech=tech.var,
                     p.value=pval, FDR=p.adjust(pval, method="BH"),
                     row.names=gnames)
    metadata(out) <- list(num.cells=ncells, resid.df=resid.df)
    return(out)
}

setGeneric("decomposeVar", function(x, fit, ...) standardGeneric("decomposeVar"))

setMethod("decomposeVar", c("ANY", "list"), .decompose_var)

setMethod("decomposeVar", c("SingleCellExperiment", "list"), 
          function(x, fit, subset.row=NULL, ..., assay.type="logcounts", get.spikes=NA) {
              
    # Checking whether we want to retrieve spikes but not use them in FDR calculations.
    if (is.na(get.spikes)) { 
        get.spikes <- TRUE
        names.to.kill <- rownames(x)[isSpike(x)]
    } else {
        names.to.kill <- character(0)
    }

    subset.row <- .SCE_subset_genes(subset.row, x=x, get.spikes=get.spikes)
    .check_centered_SF(x, assay.type=assay.type)
    out <- .decompose_var(assay(x, i=assay.type), fit, ..., subset.row=subset.row)

    # Wiping out the p-values for the spike-in transcripts, if requested. 
    if (length(names.to.kill)) {
        m <- match(names.to.kill, rownames(out))
        out$p.value[m] <- NA_real_
        out$FDR <- p.adjust(out$p.value, method="BH")
    }
    return(out)
}) 


