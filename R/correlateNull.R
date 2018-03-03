#' @export
correlateNull <- function(ncells, iters=1e6, block=NULL, design=NULL, residuals=FALSE) 
# This builds a null distribution for the modified Spearman's rho.
#
# written by Aaron Lun
# created 10 February 2016    
{
    if (!is.null(block)) { 
        if (!missing(ncells)) { 
            stop("cannot specify both 'ncells' and 'block'")
        }
        groupings <- table(block)

        # Estimating the correlation as a weighted mean of the correlations in each group.
        # This avoids the need for the normality assumption in the residual effect simulation.
        out <- 0
        for (ngr in groupings) {
            out.g <- .Call(cxx_get_null_rho, ngr, as.integer(iters))
            out <- out + out.g * ngr
        }
        out <- out/length(block)
        attrib <- list(block=block)
        
    } else if (!is.null(design)) { 
        if (!missing(ncells)) { 
            stop("cannot specify both 'ncells' and 'design'")
        }

        if (residuals) { 
            .Deprecated(msg="'residuals=TRUE' is deprecated, choose between 'design' and 'block'")
        } else {
            groupings <- .is_one_way(design)
            if (!is.null(groupings)) {
                .Deprecated(msg="'residuals=FALSE' for one-way layouts is deprecated, use 'block'")
            }
        }
        
        # Using residualsd residual effects if the design matrix is not a one-way layout (or if forced by residuals=TRUE).
        QR <- .ranksafe_qr(design)
        out <- .Call(cxx_get_null_rho_design, QR$qr, QR$qraux, as.integer(iters))
        attrib <- list(design=design)

    } else {
        out <- .Call(cxx_get_null_rho, as.integer(ncells), as.integer(iters))
        attrib <- NULL
    }

    # Storing attributes, to make sure it matches up.
    out <- sort(out)
    attributes(out) <- attrib
    return(out)  
}

