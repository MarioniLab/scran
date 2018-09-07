#' @export
#' @importFrom stats nls residuals median coef
cleanSizeFactors <- function(size.factors, num.detected, control=nls.control(warnOnly=TRUE), iterations=3, nmads=3, ...) {
    keep <- size.factors > 0
    if (all(keep)) {
        return(size.factors)
    }

    if (length(size.factors)!=length(num.detected)) {
        stop("'size.factors' and 'num.detected' should be the same length")
    }
    X <- size.factors[keep]
    Y <- num.detected[keep]
    if (length(X) < 3) {
        stop("need at least three positive values for trend fitting")
    }

    # Robustifying iterations with tricube weights.
    init <- c(A=1, B=1/max(Y))
    init <- log(init) # using log to guarantee positivity.
    weights <- rep(1, length(Y))
    lY <- log(Y)

    for (i in seq_len(iterations+1L)) {
        # Log-transforming both sides ensure positivity, avoid large libraries dominating the least-squares.
        fit <- nls(lY ~ A + log(X) - log(X * exp(B) + 1), start=init, weights=weights, control=control, ...)
        init <- coef(fit)

        resids <- abs(residuals(fit))
        bandwidth <- pmax(1e-8, median(resids) * nmads)
        weights <- (1-pmin(resids/bandwidth, 1)^3)^3
    }

    # Fitting negative size factors to this trend.
    failed <- num.detected[!keep]
    coefs <- exp(coef(fit))
    new.sf <- failed/(coefs["A"] - coefs["B"] * failed)

    # If failed > A/B (very unlikely), we set it to the maximum SF. 
    new.sf[new.sf < 0] <- max(X)

    size.factors[!keep] <- new.sf
    size.factors
}
