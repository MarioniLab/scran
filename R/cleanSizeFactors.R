#' @export
#' @importFrom stats nls residuals median coef lm
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

    # Solving for initial conditions.
    lower <- median(size.factors) > size.factors
    fit <- lm(num.detected[lower] ~ 0 + size.factors[lower])
    A <- unname(coef(fit))

    below <- which(Y < A*X)
    top <- below[which.max(Y[below])]
    B <- A / Y[top] - 1 / X[top] # must be positive due to 'below'.

    # Using logs to enforce positivity.
    init <- log(c(logA=A, logB=B)) 
    lY <- log(Y)
    lX <- log(X)

    # Robustifying iterations with tricube weights.
    weights <- rep(1, length(Y))
    for (i in seq_len(iterations+1L)) {
        # Log-transforming both sides avoids large libraries dominating the least-squares.
        fit <- nls(lY ~ logA + lX - log(X * exp(logB) + 1), start=init, weights=weights, control=control, ...)
        init <- coef(fit)

        resids <- abs(residuals(fit))
        bandwidth <- pmax(1e-8, median(resids) * nmads)
        weights <- (1-pmin(resids/bandwidth, 1)^3)^3
    }

    # Fitting negative size factors to this trend.
    failed <- num.detected[!keep]
    coefs <- coef(fit)
    new.sf <- failed/(exp(coefs["logA"]) - exp(coefs["logB"]) * failed)

    # If any failed > A/B (very unlikely), we set it to the maximum positive SF. 
    new.sf[new.sf < 0] <- max(X)

    size.factors[!keep] <- new.sf
    size.factors
}
