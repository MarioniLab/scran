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

    # Handles solo inputs.
    expect_equal(p1, combinePValues(p1, method=method, weights=weights[1]))
    expect_equal(p2, combinePValues(p2, method=method, weights=weights[2]))
    expect_equal(p3, combinePValues(p3, method=method, weights=weights[3]))

    # Handles empty inputs.
    expect_equal(
        combinePValues(p1[0], p2[0], p3[0], method=method, weights=weights),
        numeric(0))

    # Throws on invalid inputs.
    expect_error(combinePValues(p1, p2[0], method=method, weights=weights), "should have the same length")
        
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
    
    # Behaves sensibly at edge cases.
    expect_equal(combinePValues(0, 0, method="fisher"), 0)
    expect_equal(combinePValues(0, 1, method="fisher"), 0)
    expect_equal(combinePValues(1, 1, method="fisher"), 1)
})

test_that("Simes' method works correctly", {
    pout <- combinePValues(p1, p2, p3, method="simes")
    expect_equal(pout, apply(cbind(p1, p2, p3), 1, FUN=function(p) { min(p.adjust(p, method="BH")) })) 
    TESTER(p1, p2, p3, method="simes")

    # Handles ties correctly.
    pout <- combinePValues(p1, p2, p1, method="simes")
    expect_equal(pout, apply(cbind(p1, p2, p1), 1, FUN=function(p) { min(p.adjust(p, method="BH")) })) 
    TESTER(p1, p2, p1, method="simes")

    pout <- combinePValues(p1, p1, p1, method="simes")
    expect_equal(pout, apply(cbind(p1, p1, p1), 1, FUN=function(p) { min(p.adjust(p, method="BH")) })) 
    TESTER(p1, p1, p1, method="simes")

    # Behaves sensibly at edge cases.
    expect_equal(combinePValues(0, 0, method="simes"), 0)
    expect_equal(combinePValues(0, 1, method="simes"), 0)
    expect_equal(combinePValues(1, 1, method="simes"), 1)
})

test_that("Middle-Holm method works correctly", {
    # Default midway settings.
    pout <- combinePValues(p1, p2, p3, method="holm-middle")
    expect_equal(pout, apply(cbind(p1, p2, p3), 1, FUN=function(p) { median(p.adjust(p, method="holm")) })) 
    TESTER(p1, p2, p3, method="holm-middle")

    p4 <- (p1+p2)/2
    pout <- combinePValues(p1, p2, p3, p4, method="holm-middle")
    expect_equal(pout, apply(cbind(p1, p2, p3, p4), 1, FUN=function(p) { sort(p.adjust(p, method="holm"))[2] })) 
    TESTER(p1, p2, p3, method="holm-middle")

    # Testing alternative min.prop settings.
    pout <- combinePValues(p1, p2, p3, p4, method="holm-middle", min.prop=0.2)
    expect_equal(pout, apply(cbind(p1, p2, p3, p4), 1, FUN=function(p) { sort(p.adjust(p, method="holm"))[1] })) 

    pout <- combinePValues(p1, p2, p3, p4, method="holm-middle", min.prop=0.25)
    expect_equal(pout, apply(cbind(p1, p2, p3, p4), 1, FUN=function(p) { sort(p.adjust(p, method="holm"))[1] })) 
    
    pout <- combinePValues(p1, p2, p3, p4, method="holm-middle", min.prop=0.45)
    expect_equal(pout, apply(cbind(p1, p2, p3, p4), 1, FUN=function(p) { sort(p.adjust(p, method="holm"))[2] })) 

    pout <- combinePValues(p1, p2, p3, p4, method="holm-middle", min.prop=0.65)
    expect_equal(pout, apply(cbind(p1, p2, p3, p4), 1, FUN=function(p) { sort(p.adjust(p, method="holm"))[3] })) 

    pout <- combinePValues(p1, p2, p3, p4, method="holm-middle", min.prop=0.75)
    expect_equal(pout, apply(cbind(p1, p2, p3, p4), 1, FUN=function(p) { sort(p.adjust(p, method="holm"))[3] })) 

    pout <- combinePValues(p1, p2, p3, p4, method="holm-middle", min.prop=0.85)
    expect_equal(pout, apply(cbind(p1, p2, p3, p4), 1, FUN=function(p) { sort(p.adjust(p, method="holm"))[4] })) 

    # Handles ties correctly.
    pout <- combinePValues(p1, p2, p1, method="holm-middle")
    expect_equal(pout, apply(cbind(p1, p2, p1), 1, FUN=function(p) { median(p.adjust(p, method="holm")) })) 
    TESTER(p1, p2, p1, method="simes")

    pout <- combinePValues(p1, p1, p1, method="holm-middle")
    expect_equal(pout, apply(cbind(p1, p1, p1), 1, FUN=function(p) { median(p.adjust(p, method="holm")) })) 
    TESTER(p1, p1, p1, method="simes")

    # Behaves sensibly at edge cases.
    expect_equal(combinePValues(0, 0, method="holm-middle"), 0)
    expect_equal(combinePValues(0, 1, method="holm-middle"), 0)
    expect_equal(combinePValues(1, 1, method="holm-middle"), 1)
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

    # Handles weights as a list.
    expect_equal(pnorm(Q), combinePValues(p1, p2, p3, method="z", weights=lapply(W, rep, length.out=length(p1))))

    # Throws errors correctly.
    expect_error(combinePValues(p1,p2,p3,method="z", weights=1), "should be equal")
    expect_error(combinePValues(p1,p2,p3,method="z", weights=c(-1, 1, 2)), "must be positive")
    expect_error(combinePValues(p1,p2,p3,method="z", weights=c(NA, 1, 2)), "must be positive")

    # Behaves sensibly at edge cases.
    expect_equal(combinePValues(0, 0, method="z"), 0)
    expect_equal(combinePValues(0, 1, method="z"), 0.5)
    expect_equal(combinePValues(1, 1, method="z"), 1)
})

test_that("Berger's IUT works correctly", {
    expect_equal(combinePValues(p1, p2, p3, method="berger"), pmax(p1, p2, p3)) 
    TESTER(p1, p2, p3, method="berger")

    # Behaves sensibly at edge cases.
    expect_equal(combinePValues(0, 0, method="berger"), 0)
    expect_equal(combinePValues(0, 1, method="berger"), 1)
    expect_equal(combinePValues(1, 1, method="berger"), 1)
})
