.decompose_var <- function(x, fit, design=NA, subset.row=NULL, ...)
# Computes the biological variability of the log-CPMs by subtracting the
# inferred technical variance from the total variance.
#
# written by Aaron Lun
# created 21 January 2016 
# last modified 11 July 2017
{
    if (!missing(x)) { 
        subset.row <- .subset_to_index(subset.row, x, byrow=TRUE)
        checked <- .make_var_defaults(x, fit=fit, design=design)
        design <- checked$design
        QR <- .ranksafe_qr(design)
        
        lout <- .Call(cxx_fit_linear_model, QR$qr, QR$qraux, x, subset.row-1L, FALSE)
        lmeans <- lout[[1]]
        lvar <- lout[[2]]
        gnames <- rownames(x)[subset.row]
    } else {
        lmeans <- fit$mean
        lvar <- fit$var
        design <- fit$design
        gnames <- names(lmeans)
    }
    
    tech.var <- fit$trend(lmeans)
    bio.var <- lvar - tech.var
    pval <- testVar(total=lvar, null=tech.var, df=nrow(design) - ncol(design), second.df=fit$df2, ...)
    out <- data.frame(mean=lmeans, total=lvar, bio=bio.var, tech=tech.var,
                      p.value=pval, FDR=p.adjust(pval, method="BH"),
                      row.names=gnames)
    return(out)
}

setGeneric("decomposeVar", function(x, fit, ...) standardGeneric("decomposeVar"))

setMethod("decomposeVar", c("ANY", "list"), .decompose_var)

setMethod("decomposeVar", c("SingleCellExperiment", "list"), 
          function(x, fit, subset.row=NULL, ..., assay.type="logcounts", get.spikes=FALSE) {

    subset.row <- .SCE_subset_genes(subset.row, x=x, get.spikes=get.spikes)
    .check_centered_SF(x, assay.type=assay.type)
    .decompose_var(assay(x, i=assay.type), fit, ..., subset.row=subset.row)
})

