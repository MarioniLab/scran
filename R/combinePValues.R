#' @export
#' @importFrom stats pnorm qnorm pchisq
combinePValues <- function(..., method=c("fisher", "z", "simes", "berger"), weights=NULL, log.p=FALSE) 
# Combines p-values with a variety of methods.
# Weights are support in Stouffer's Z method.
#
# written by Aaron Lun
# created 16 September 2018
{
    input <- list(...)
    Np <- unique(lengths(input))
    if (length(Np) != 1) {
        stop("all p-value vectors must have the same length")
    }

    method <- match.arg(method)
    switch(method,
        fisher={
            if (log.p) {
                all.logp <- input
            } else {
                all.logp <- lapply(input, FUN=log)
            }

            n <- integer(Np)
            X <- numeric(Np)
            for (i in seq_along(all.logp)) {
                current <- all.logp[[i]]
                keep <- !is.na(current)
                X[keep] <- X[keep] + current[keep]
                n <- n + keep
            }
            
            n[n==0] <- NA_real_ # ensure that we get NA outputs.
            pchisq(-2*X, df=2*n, lower.tail=FALSE, log.p=log.p)
        },
        simes={
            .Call(cxx_combine_simes, input, log.p)
        }, 
        z={
            if (is.null(weights)) {
                weights <- rep(1, length(input))
            } else if (length(weights)!=length(input)) {
                stop("'length(weights)' must be equal to number of vectors in '...'")
            } else {
                check <- unlist(lapply(weights, range))
                if (any(is.na(check) | check<=0)) {
                    stop("weights must be positive")
                }
            }

            Z <- W2 <- numeric(Np)
            for (i in seq_along(input)) {
                current <- input[[i]]
                keep <- !is.na(current)
                Z[keep] <- Z[keep] + qnorm(current[keep], log.p=log.p) * weights[[i]]
                W2 <- W2 + ifelse(keep, weights[[i]]^2, 0)
            }

            # Combining p-values of 0 and 1 will yield zscores of -Inf + Inf => NaN.
            # Here, we set them to 0 to get p-values of 0.5, as the Z-method doesn't
            # give coherent answers when you have p-values at contradicting extremes.
            Z[is.nan(Z)] <- 0

            W2[W2==0] <- NA_real_ # ensure we get NA outputs.
            pnorm(Z/sqrt(W2), log.p=log.p) 
        },
        berger={
            do.call(pmax, c(input, list(na.rm=TRUE)))
        }
    )
}
