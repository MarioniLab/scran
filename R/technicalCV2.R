.technicalCV2 <- function(x, is.spike, sf.cell=NULL, sf.spike=NULL, 
                          cv2.limit=0.3, cv2.tol=0.8, min.bio.disp=0.25) 
# This implements Brennecke's method for fitting spike-ins to the CV2 of various genes.
# The output can also be plotted fairly easily.
#
# written by Aaron Lun
# based on code by Phillipe Brennecke et al. (2013).
# created 11 July 2016
{
    if (any(!is.na(is.spike))) { 
        if (any(is.na(is.spike))) { 
            stop("missing values in 'is.spike'")
        }
        is.spike <- .subset_to_index(is.spike, x, byrow=TRUE)
        is.cell <- seq_len(nrow(x))[-is.spike]
    } else {
        is.cell <- is.spike <- seq_len(nrow(x))
    }
    if (length(is.spike) < 2L) {
        stop("need at least 2 spike-ins for trend fitting")
    }

    # Computing size factors, if not supplied.
    if (is.null(sf.cell)) { 
        if (length(is.cell)) { 
            sf.cell <- DESeq2::estimateSizeFactorsForMatrix(x[is.cell,,drop=FALSE])
        } else {
            sf.cell <- rep(1, ncol(x)) # Any value will do here.
        }
    } 
    if (is.null(sf.spike)) {
        sf.spike <- DESeq2::estimateSizeFactorsForMatrix(x[is.spike,,drop=FALSE])
    }

    # Computing the statistics.
    cell.out <- .Call(cxx_compute_CV2, x, is.cell-1L, sf.cell)
    if (is.character(cell.out)) { 
        stop(cell.out)
    }
    spike.out <- .Call(cxx_compute_CV2, x, is.spike-1L, sf.spike)
    if (is.character(spike.out)) {
        stop(spike.out)
    }

    means.cell <- cell.out[[1]]
    vars.cell <- cell.out[[2]]
    cv2.cell <- vars.cell/means.cell^2

    means.spike <- spike.out[[1]]
    vars.spike <- spike.out[[2]]
    cv2.spike <- vars.spike/means.spike^2

    # Fitting the trend.
    above.limit <- cv2.spike > cv2.limit
    if (any(above.limit)) { 
        minMeanForFitA <- unname( quantile( means.spike[ above.limit ], cv2.tol ) )
        useForFitA <- means.spike >= minMeanForFitA
    } else {
        warning("no spike-ins above 'cv2.limit', using all spikes for trend fitting")
        useForFitA <- !logical(length(means.spike)) # using all spikes.
    }
    fitA <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/means.spike[useForFitA] ), cv2.spike[useForFitA] )

    #  Testing whether the null is true.
    xi <- mean( 1 / sf.spike )
    m <- ncol(x)
    psia1thetaA <- mean( 1 / sf.spike ) + ( coefficients(fitA)["a1tilde"] - xi ) * mean( sf.spike / sf.cell )
    cv2thA <- coefficients(fitA)["a0"] + min.bio.disp + coefficients(fitA)["a0"] * min.bio.disp
    testDenomA <- ( means.cell * psia1thetaA + means.cell^2 * cv2thA ) / ( 1 + cv2thA/m )
    pA <- 1 - pchisq( vars.cell * (m-1) / testDenomA, m-1 )

    # Formatting the returned output.
    output.mean <- output.var <- output.cv2 <- numeric(nrow(x))
    output.mean[is.spike] <- means.spike
    output.var[is.spike] <- vars.spike
    output.cv2[is.spike] <- cv2.spike
    output.mean[is.cell] <- means.cell
    output.var[is.cell] <- vars.cell
    output.cv2[is.cell] <- cv2.cell
    output.trend <- coefficients(fitA)["a0"] + coefficients(fitA)["a1tilde"]/output.mean
    output.p <- rep(NA_real_, nrow(x))
    output.p[is.cell] <- pA
    return(data.frame(mean=output.mean, var=output.var, cv2=output.cv2, trend=output.trend,
                      p.value=output.p, FDR=p.adjust(output.p, method="BH"), row.names=rownames(x)))
}

setGeneric("technicalCV2", function(x, ...) standardGeneric("technicalCV2"))

setMethod("technicalCV2", "matrix", .technicalCV2)

setMethod("technicalCV2", "SCESet", function(x, spike.type=NULL, ..., assay="counts") {
    sf.cell <- sizeFactors(x)

    if (is.null(spike.type) || !is.na(spike.type)) {  
        is.spike <- isSpike(x, type=spike.type)
        if (is.null(spike.type)) { 
            # Get all spikes.
            spike.type <- whichSpike(x)            
        }
        if (!length(spike.type)) { 
            stop("no spike-in sets specified from 'x'")
        }

        # Collecting the size factors for the requested spike-in sets.
        collected <- list() 
        for (st in spike.type) {
            cur.sf <- suppressWarnings(sizeFactors(x, type=st))
            if (is.null(cur.sf)) { 
                collected[st] <- list(NULL)
            } else {
                collected[[st]] <- cur.sf
            }
        }

        # Check that all spike-in factors are either NULL or identical.
        if (length(collected)) {
            for (i in seq_along(collected)[-1]) { 
                if (!isTRUE(all.equal(collected[[i]], collected[[1]]))) { 
                    stop("size factors differ between spike-in sets")
                }
            }
            sf.spike <- collected[[1]]
        } 

        # Diverting to the cell-based size factors if all spike-in factors are NULL.
        if (is.null(sf.spike)) {
            warning("no spike-in size factors set, using cell-based factors")
            sf.spike <- sf.cell
        }

    } else {
        sf.spike <- sf.cell
        is.spike <- NA
    }
    .technicalCV2(assayDataElement(x, assay), is.spike=is.spike, 
                  sf.cell=sf.cell, sf.spike=sf.spike, ...)          
})

