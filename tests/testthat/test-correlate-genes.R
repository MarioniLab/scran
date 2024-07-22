# Tests the correlateGenes function.
# library(testthat); library(scran); source("test-correlate-genes.R")

set.seed(100022)
Ngenes <- 20
Ncells <- 100
X <- log(matrix(rpois(Ngenes*Ncells, lambda=10), nrow=Ngenes) + 1)
rownames(X) <- paste0("X", seq_len(Ngenes))
ref <- correlatePairs(X)

test_that("correlateGenes works correctly", {
    out <- correlateGenes(ref)
    for (x in rownames(X)) {
        collected <- ref$gene1 == x | ref$gene2==x

        simes.p <- min(p.adjust(ref$p.value[collected], method="BH"))
        expect_equal(simes.p, out$p.value[out$gene==x])

        max.i <- which.max(abs(ref$rho[collected]))
        expect_equal(ref$rho[collected][max.i], out$rho[out$gene==x])
    }
})

test_that("correlateGenes handles silly inputs", {
    out <- correlateGenes(ref[1:10,])
    expect_identical(out[0,], correlateGenes(ref[0,]))
})
