# This tests the testLinearModel function.
# library(testthat); library(scran); source("test-linear-test.R")

test_that("linear model testing works for contrast vectors in categorical designs", {
    y <- matrix(rnorm(10000), ncol=100)
                                                 
    A <- gl(2, 50)
    design <- model.matrix(~A)
    out <- testLinearModel(y, design, contrast=c(0, 1))

    for (i in seq_len(10)) {
        Y <- y[i,]
        fit <- lm(Y ~ A)
        stats <- summary(fit)$coefficients
        expect_equal(stats[2,"Pr(>|t|)"], out$p.value[i])
    }

    # Robust to another formulation of the same comparison:
    design2 <- model.matrix(~0 + A)
    alt <- testLinearModel(y, design2, contrast=c(-1, 1))
    expect_equal(out, alt)

    alt <- testLinearModel(y, design)
    expect_equal(out, alt)
})

test_that("linear model testing works for contrast vectors with continuous variables", {
    y <- matrix(rnorm(10000), ncol=100)

    u <- runif(100)
    design <- model.matrix(~u)
    out <- testLinearModel(y, design, contrast=c(0, 1))

    for (i in seq_len(10)) {
        Y <- y[i,]
        fit <- lm(Y ~ u)
        stats <- summary(fit)$coefficients
        expect_equal(stats[2,"Pr(>|t|)"], out$p.value[i])
    }

    alt <- testLinearModel(y, design)
    expect_equal(out, alt)

    # Handles more complex designs like a champ.
    v <- gl(4, 25)
    design <- model.matrix(~u + v)
    out <- testLinearModel(y, design, contrast=c(0, 1, 0, 0, 0))

    for (i in seq_len(10)) {
        Y <- y[i,]
        fit <- lm(Y ~ u + v)
        stats <- summary(fit)$coefficients
        expect_equal(stats[2,"Pr(>|t|)"], out$p.value[i])
    }

    alt <- testLinearModel(y, design, coef=2)
    expect_equal(out, alt)
})

test_that("linear model testing works for contrast matrices", {
    y <- matrix(rnorm(10000), ncol=100)

    A <- gl(2, 50)
    u <- runif(100)
    design <- model.matrix(~A + u)
    out <- testLinearModel(y, design, contrast=cbind(c(0, 1, 0), c(0, 0, 1)))

    for (i in seq_len(10)) {
        Y <- y[i,]
        fit <- lm(Y ~ A + u)
        fit0 <- lm(Y ~ 1)
        stats <- anova(fit, fit0)
        expect_equal(stats[2,"Pr(>F)"], out$p.value[i])
    }

    alt <- testLinearModel(y, design, coef=2:3)
    expect_equal(out, alt)

    # Another design:
    B <- gl(4, 25)
    design <- model.matrix(~B)
    out <- testLinearModel(y, design, contrast=cbind(c(0,1,0,0), c(0,0,1,-1)))

    for (i in seq_len(10)) {
        Y <- y[i,]
        fit <- lm(Y ~ B)
        fit0 <- lm(Y ~ A)
        stats <- anova(fit, fit0)
        expect_equal(stats[2,"Pr(>F)"], out$p.value[i])
    }
})

test_that("linear model testing works with miscellaneous options", {
    y <- matrix(rnorm(10000), ncol=100)

    A <- gl(2, 50)
    u <- runif(100)
    design <- model.matrix(~A + u)
    ref <- testLinearModel(y, design, contrast=cbind(c(0, 1, 0), c(0, 0, 1)))
    out <- testLinearModel(y, design, contrast=cbind(c(0, 1, 0), c(0, 0, 1)), subset.row=1:10)

    ref$FDR <- out$FDR <- NULL
    expect_identical(ref[1:10,], out)

    # Respects row names.
    rownames(y) <- paste0("GENE_", seq_len(nrow(y)))
    out <- testLinearModel(y, design, contrast=cbind(c(0, 1, 0), c(0, 0, 1)))
    expect_identical(rownames(out), rownames(y))

    out <- testLinearModel(y, design, contrast=cbind(c(0, 1, 0), c(0, 0, 1)), subset.row=1:10)
    expect_identical(rownames(out), rownames(y)[1:10])
})

