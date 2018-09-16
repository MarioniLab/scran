# This tests the combinePValues function.
# library(scran); library(testthat); source("test-combine-p.R")

set.seed(20000)
p1 <- runif(1000)
p2 <- runif(1000)
p3 <- runif(1000)    

TESTER <- function(p1, p2, p3, method, weights=NULL) {
    # Behaves properly upon logging.
    expect_equal(
        combinePValues(log(p1), log(p2), log(p3), method=method, weights=weights, log.p=TRUE),
        log(combinePValues(p1, p2, p3, method=method, weights=weights))
    )

    # Handles empty inputs.
    expect_equal(
        combinePValues(p1[0], p2[0], p3[0], method=method, weights=weights),
        numeric(0))

    # Throws on invalid inputs.
    expect_error(combinePValues(p1, p2[0], method=method, weights=weights), "must have the same length")
        
    # Handles partial NA values correctly.
    some.na <- sample(length(p1), length(p1)/2)
    p1.na <- p1
    p1.na[some.na] <- NA
    p2.na <- p2
    p2.na[-some.na] <- NA

    ref <- combinePValues(p1, p3, method=method, weights=weights[c(1,3)])
    ref2 <- combinePValues(p2, p3, method=method, weights=weights[c(2,3)])
    ref[some.na] <- ref2[some.na]
    expect_equal(ref, combinePValues(p1.na, p2.na, p3, method=method, weights=weights))

    # Handles all-NA values correctly.
    p1.na <- p1
    p2.na <- p2
    p3.na <- p3
    p1.na[1] <- p2.na[1] <- p3.na[1] <- NA_real_
    out <- combinePValues(p1.na, p2.na, p3.na, method=method, weights=weights)
    expect_equal(out[1], NA_real_)
    expect_false(any(is.na(out[-1])))

    return(TRUE)
}

test_that("Fisher's method works correctly", {
    pout <- combinePValues(p1, p2, p3, method="fisher")
    expect_equal(pout, pchisq(-2*rowSums(log(cbind(p1, p2, p3))), df=6, lower.tail=FALSE))
    TESTER(p1, p2, p3, method="fisher")
})

test_that("Simes' method works correctly", {
    pout <- combinePValues(p1, p2, p3, method="simes")
    expect_equal(pout, apply(cbind(p1, p2, p3), 1, FUN=function(p) { min(p.adjust(p, method="BH")) })) 
    TESTER(p1, p2, p3, method="simes")
})

test_that("Stouffer's Z method works correctly", {
    Q <- qnorm(rbind(p1, p2, p3))
    Q <- colSums(Q)/sqrt(nrow(Q))
    pout <- combinePValues(p1, p2, p3, method="z")
    expect_equal(pout, pnorm(Q))
    
    TESTER(p1, p2, p3, method="z")

    # Handles weights correctly.
    expect_equal(pout, combinePValues(p1, p2, p3, method="z", weights=c(2,2,2)))
    
    Q <- qnorm(rbind(p1, p2, p3))
    W <- 1:3
    Q <- colSums(Q * W)/sqrt(sum(W^2))
    expect_equal(pnorm(Q), combinePValues(p1, p2, p3, method="z", weights=W))

    TESTER(p1, p2, p3, method="z", weights=c(5, 2, 3))

    # Throws errors correctly.
    expect_error(combinePValues(p1,p2,p3,method="z", weights=1), "must be equal")
})

test_that("Berger's IUT works correctly", {
    expect_equal(combinePValues(p1, p2, p3, method="berger"), pmax(p1, p2, p3)) 
    TESTER(p1, p2, p3, method="berger")
})
