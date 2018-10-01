# This tests the testVar function.
# require(scran); require(testthat); source("test-testvar.R")

set.seed(20002)
true.p <- runif(100)
trended <- runif(100, 1, 2)
df <- 20
observed <- trended * qchisq(true.p, df=df, lower.tail=FALSE)/df

test_that("testVar's chi-squared test works as expected", {
    pvals <- testVar(observed, trended, df=df)
    expect_equal(pvals, true.p)
    pvals <- testVar(observed, trended, df=df, log.p=TRUE)
    expect_equal(pvals, log(true.p))

    design <- model.matrix(~factor(rep(c(1,2), each=11)))
    pvals <- testVar(observed, trended, design=design)
    expect_equal(pvals, true.p)
})

test_that("testVar's F-test works as expected", {
    nspikes <- ncells <- 100
    spike.means <- 2^runif(nspikes, -2, 8)
    spike.disp <- (100/spike.means + 0.5) * 1/rchisq(nspikes, df=10)
    spike.data <- matrix(rnbinom(nspikes*ncells, mu=spike.means, size=1/spike.disp), ncol=ncells)

    exprs <- log2(spike.data/(colSums(spike.data)/mean(colSums(spike.data)))+1)
    fit <- trendVar(exprs)
    pvals <- testVar(fit$var, fit$trend(fit$mean), df=ncells-1, second.df=fit$df2, test='f')

    rat <- (fit$var/fit$trend(fit$mean))
    df1 <- ncells - 1L
    ffit <- limma::fitFDistRobustly(rat[fit$var > 0 & fit$mean >= 0.1], df=df1) # filtering out zero-variance, low-abundance genes.
    expect_equal(pvals, pf(rat/ffit$scale, df1=df1, df2=ffit$df2, lower.tail=FALSE))

    # Testing log-transformation.
    lpvals <- testVar(fit$var, fit$trend(fit$mean), df=ncells-1, second.df=fit$df2, test='f', log.p=TRUE)
    expect_equal(log(pvals), lpvals)

    # Checking for error upon no df.2
    expect_error(testVar(fit$var, fit$trend(fit$mean), df=ncells-1, test='f'),
                 "second df from trendVar() must be specified for test='f'", fixed=TRUE)
})

test_that("testVar works with silly inputs", {
    expect_identical(testVar(0, 0, df=10), 1)
    expect_identical(testVar(0, 0, df=10, test="f", second.df=2), 1)
    expect_identical(testVar(0, 0, df=0), NA_real_)
    expect_identical(testVar(0, 0, df=0, test="f", second.df=2), NA_real_)

    expect_identical(testVar(numeric(0), trended, df=df), numeric(0))
    expect_identical(testVar(numeric(0), trended, df=df, test="f", second.df=2), numeric(0))

    expect_identical(testVar(observed, numeric(0), df=df), numeric(0))
    expect_identical(testVar(observed, numeric(0), df=df, test="f", second.df=2), numeric(0))

    expect_identical(testVar(observed, trended, df=numeric(0)), rep(NA_real_, length(observed)))
    expect_identical(testVar(observed, trended, df=numeric(0), test="f", second.df=2), rep(NA_real_, length(observed)))
})
