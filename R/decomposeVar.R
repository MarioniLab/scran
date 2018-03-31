#' @importFrom S4Vectors DataFrame "metadata<-"
#' @importFrom stats qnorm pnorm p.adjust
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
        # Setting up the input for combineVar.
        nbatches <- ncol(lmeans)
        collected <- vector("list", nbatches) 
        names(collected) <- colnames(lmeans)

        for (x in names(collected)) { 
            collected[[x]] <- .create_var_df(lmeans[,x], lvar[,x], fit$trend, resid.df=resid.df[,x], 
                                             gene.names=rownames(lmeans), ncells=sum(block==x), second.df=fit$df2)
        }
        
        out <- do.call(combineVar, c(collected, list(method="z")))
        return(out)
    } 
        
    .create_var_df(lmeans, lvar, fit$trend, resid.df=resid.df, 
                   gene.names=names(lmeans), ncells=ncells, second.df=fit$df2)
}

.create_var_df <- function(lmean, lvar, trendfun, resid.df, gene.names, ncells, ...) 
# Helper function to create the output DataFrame. 
{
    lvar <- as.numeric(lvar)
    tech.var <- as.numeric(trendfun(lmean))
    pval <- testVar(total=lvar, null=tech.var, df=unname(resid.df), ...)
    bio.var <- lvar - tech.var
    out <- DataFrame(mean=as.numeric(lmean), 
                     total=lvar, bio=bio.var, tech=tech.var,
                     p.value=pval, FDR=p.adjust(pval, method="BH"),
                     row.names=gene.names)
    metadata(out) <- list(num.cells=ncells, resid.df=resid.df)
    return(out)
}

##############################
# Defining the S4 internals. #
##############################

#' @export
setGeneric("decomposeVar", function(x, fit, ...) standardGeneric("decomposeVar"))

#' @export
setMethod("decomposeVar", c("ANY", "list"), .decompose_var)

#' @importFrom SummarizedExperiment assay
#' @importFrom stats p.adjust
#' @export
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


