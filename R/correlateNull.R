correlateNull <- function(ncells, iters=1e6, design=NULL, residuals=FALSE) 
# This builds a null distribution for the modified Spearman's rho.
#
# written by Aaron Lun
# created 10 February 2016
# last modified 19 January 2016
{
    if (!is.null(design)) { 
        if (!missing(ncells)) { 
            stop("cannot specify both 'ncells' and 'design'")
        }

        groupings <- .isOneWay(design)
        if (is.null(groupings) || residuals) { 
            # Using residualsd residual effects if the design matrix is not a one-way layout (or if forced by residuals=TRUE).
            QR <- qr(design, LAPACK=TRUE)
            out <- .Call(cxx_get_null_rho_design, QR$qr, QR$qraux, as.integer(iters))
            if (is.character(out)) { 
                stop(out)
            }
        } else {
            # Otherwise, estimating the correlation as a weighted mean of the correlations in each group.
            # This avoids the need for the normality assumption in the residual effect simulation.
            out <- 0
            for (gr in groupings) {
                out.g <- .Call(cxx_get_null_rho, length(gr), as.integer(iters))
                if (is.character(out.g)) { 
                    stop(out.g)
                }
                out <- out + out.g * length(gr)
            }
            out <- out/nrow(design)
        }
        attrib <- list(design=design, residuals=residuals)

    } else {
        out <- .Call(cxx_get_null_rho, as.integer(ncells), as.integer(iters))
        if (is.character(out)) { 
            stop(out)
        }
        attrib <- NULL
    }

    # Storing attributes, to make sure it matches up.
    out <- sort(out)
    attributes(out) <- attrib
    return(out)  
}

