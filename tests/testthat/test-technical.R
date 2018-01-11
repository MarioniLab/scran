# Checks the technicalCV2 and improvedCV2 functions.
# require(scran); require(testthat); source("test-technical.R")

set.seed(6000)
ngenes <- 10000
means <- 2^runif(ngenes, 6, 10)
dispersions <- 10/means + 0.2
nsamples <- 50

counts <- matrix(rnbinom(ngenes*nsamples, mu=means, size=1/dispersions), ncol=nsamples)
sf <- 2^rnorm(nsamples)
is.spike <- logical(ngenes)
is.spike[seq_len(500)] <- TRUE

test_that("CV2 calculations are correct", {
    chosen <- which(is.spike)
    stuff <- .Call(scran:::cxx_compute_CV2, counts, chosen - 1L, sf, NULL)
    normed <- t(t(counts[chosen,])/sf)
    expect_equal(stuff[[1]], rowMeans(normed))
    expect_equal(stuff[[2]], apply(normed, 1, var))
    
    chosen <- sample(ngenes, 1000)
    stuff <- .Call(scran:::cxx_compute_CV2, counts, chosen - 1L, sf, NULL)
    normed <- t(t(counts[chosen,])/sf)
    expect_equal(stuff[[1]], rowMeans(normed))
    expect_equal(stuff[[2]], apply(normed, 1, var))

    # Testing its ability to antilog.
    prior.count <- 1
    logged <- log2(t(t(counts)/sf) + prior.count) 
    restuff <- .Call(scran:::cxx_compute_CV2, logged, chosen - 1L, NULL, prior.count)
    expect_equal(stuff, restuff)
})

#######################################################

# Setting up the object here for convenience.
rownames(counts) <- paste0("X", seq_len(ngenes))
colnames(counts) <- paste0("Y", seq_len(nsamples))
X <- SingleCellExperiment(list(counts=counts))
isSpike(X, "Spikes") <- is.spike

sizeFactors(X) <- sf
spike.sf <- 2^rnorm(nsamples)
sizeFactors(X, type="Spikes") <- spike.sf

test_that("technicalCV2 works sensibly with SingleCellExperiment", {
    default <- technicalCV2(counts, is.spike, sf.cell=sf, sf.spike=spike.sf)
    as.sceset <- technicalCV2(X)
    expect_equal(default, as.sceset)
    as.sceset <- technicalCV2(X, spike.type="Spikes")
    expect_equal(default, as.sceset)

    # Checking what happens with multiple spike-in sets.
    X2 <- SingleCellExperiment(list(counts=counts))
    expect_error(suppressWarnings(technicalCV2(X2)), "no spike-in sets specified from 'x'")
    subset <- split(which(is.spike), rep(1:2, length.out=sum(is.spike)))
    isSpike(X2, "MySpike") <- subset[[1]]
    isSpike(X2, "SecondSpike") <- subset[[2]]

    sizeFactors(X2) <- sf
    sizeFactors(X2, type="MySpike") <- spike.sf
    sizeFactors(X2, type="SecondSpike") <- spike.sf
    
    expect_equal(technicalCV2(X2), as.sceset)
    expect_equal(technicalCV2(X2, spike.type=c("MySpike", "SecondSpike")), as.sceset)
    expect_equal(technicalCV2(X2, spike.type="MySpike"), technicalCV2(counts, is.spike=subset[[1]], sf.cell=sf, sf.spike=spike.sf))
   
    # Checking what happens when the spike-in size factors aren't set. 
    sizeFactors(X2, type="SecondSpike") <- rev(spike.sf)
    expect_error(technicalCV2(X2), "size factors differ between spike-in sets")
    expect_equal(technicalCV2(X2, spike.type="SecondSpike"), technicalCV2(counts, is.spike=subset[[2]], sf.cell=sf, sf.spike=rev(spike.sf)))
    
    sizeFactors(X2, type="SecondSpike") <- NULL
    expect_error(technicalCV2(X2), "size factors differ between spike-in sets")
    
    sizeFactors(X2, type="MySpike") <- NULL
    expect_warning(expect_equal(technicalCV2(X2), technicalCV2(counts, is.spike=is.spike, sf.cell=sf, sf.spike=sf)),
                   "no spike-in size factors set, using cell-based factors")
})

test_that("technicalCV2 fits to all genes properly", {
    all.used <- technicalCV2(counts, is.spike=NA)
    counts2 <- rbind(counts, counts)
    rownames(counts2) <- paste0("X", seq_len(nrow(counts2)))
    is.spike2 <- rep(c(FALSE, TRUE), each=nrow(counts))
    all.used2 <- technicalCV2(counts2, is.spike2)
    expect_equal(all.used, all.used2[!is.spike2,])
   
    # Repeating for the SingleCellExperiment.
    all.used <- technicalCV2(counts(X), sf.cell=sizeFactors(X), sf.spike=sizeFactors(X), is.spike=NA)
    all.used3 <- technicalCV2(X, spike.type=NA)
    expect_equal(all.used, all.used3)
})

test_that("technicalCV2 behaves in the presence of silly inputs", {
    expect_error(technicalCV2(X, spike.type="whee"), "spike-in set 'whee' does not exist")
    expect_error(technicalCV2(X[0,], spike.type="Spikes"), "need at least 2 spike-ins for trend fitting")
    expect_error(technicalCV2(X[,0], spike.type="Spikes"), "need two or more cells to compute variances")

    out <- technicalCV2(counts, is.spike=!logical(ngenes))
    expect_true(all(is.na(out$p.value)))
    expect_error(technicalCV2(counts, is.spike=logical(ngenes)), "need at least 2 spike-ins for trend fitting")
})

# Testing for robustness to zeroes.
test_that("technicalCV2 is robust to zeroes", {
    counts[1,] <- 0
    counts[501,] <- 0
    out <- technicalCV2(counts, is.spike=NA)
    expect_identical(out$mean[1], 0)
    expect_identical(out$var[1], 0)
    expect_identical(out$cv2[1], NaN)
    expect_identical(out$trend[1], NA_real_)
    expect_equivalent(out[1,], out[501,])
})

#######################################################

X <- normalize(X)
test_that("improvedCV2 works sensibly with SingleCellExperiment", {
    default <- improvedCV2(counts, is.spike, sf.cell=sf, sf.spike=spike.sf)
    as.sceset <- improvedCV2(X)
    expect_equal(default, as.sceset)
    as.sceset <- improvedCV2(X, assay.type="counts")
    expect_equal(default, as.sceset)
    as.sceset <- improvedCV2(X, spike.type="Spikes")
    expect_equal(default, as.sceset)

    # Checking what happens with multiple spike-in sets.
    X2 <- SingleCellExperiment(list(counts=counts))
    expect_error(suppressWarnings(improvedCV2(X2)), "no spike-in sets specified from 'x'")

    subset <- split(which(is.spike), rep(1:2, length.out=sum(is.spike)))
    isSpike(X2, "MySpike") <- subset[[1]]
    isSpike(X2, "SecondSpike") <- subset[[2]]

    sizeFactors(X2) <- sf
    sizeFactors(X2, type="MySpike") <- spike.sf
    sizeFactors(X2, type="SecondSpike") <- spike.sf
    
    X2 <- normalize(X2)
    expect_equal(improvedCV2(X2), as.sceset)
    expect_equal(improvedCV2(X2, assay.type="counts"), as.sceset)
    expect_equal(improvedCV2(X2, spike.type=c("MySpike", "SecondSpike")), as.sceset)
    expect_equal(improvedCV2(X2, assay.type="counts", spike.type="MySpike"), 
                 improvedCV2(counts, is.spike=subset[[1]], sf.cell=sf, sf.spike=spike.sf))

    # Checking what happens when the spike-in size factors aren't set. 
    sizeFactors(X2, type="SecondSpike") <- rev(spike.sf)
    expect_error(improvedCV2(X2), "size factors differ between spike-in sets")
    expect_equal(improvedCV2(X2, assay.type="counts", spike.type="SecondSpike"), 
                 improvedCV2(counts, is.spike=subset[[2]], sf.cell=sf, sf.spike=rev(spike.sf)))

    sizeFactors(X2, type="SecondSpike") <- NULL
    expect_error(improvedCV2(X2, assay.type="counts"), "size factors differ between spike-in sets")

    sizeFactors(X2, type="MySpike") <- NULL
    expect_warning(expect_equal(improvedCV2(X2, assay.type="counts"), improvedCV2(counts, is.spike=is.spike, sf.cell=sf, sf.spike=sf)),
                   "no spike-in size factors set, using cell-based factors")
})

test_that("improvedCV2 fits to all genes properly", {
    all.used <- improvedCV2(counts, is.spike=NA)
    counts2 <- rbind(counts, counts)
    rownames(counts2) <- paste0("X", seq_len(nrow(counts2)))
    is.spike2 <- rep(c(FALSE, TRUE), each=nrow(counts))
    all.used2 <- improvedCV2(counts2, is.spike2)
    expect_equal(all.used, all.used2[!is.spike2,])
   
    # Repeating for the SingleCellExperiment.
    all.used <- improvedCV2(logcounts(X), log.prior=1, is.spike=NA)
    all.used3 <- improvedCV2(X, spike.type=NA)
    expect_equal(all.used, all.used3)

    all.used <- improvedCV2(counts(X), sf.cell=sizeFactors(X), sf.spike=sizeFactors(X), is.spike=NA)
    all.used3 <- improvedCV2(X, assay.type="counts", spike.type=NA)
    expect_equal(all.used, all.used3)

    # Robustness works.
    all.used <- improvedCV2(X, robust=TRUE)
})

test_that("improvedCV2 behaves in the presence of silly inputs", {
    expect_error(improvedCV2(X, spike.type="whee"), "spike-in set 'whee' does not exist")
    expect_error(improvedCV2(X[0,], spike.type="Spikes"), "need at least 2 spike-ins for trend fitting")
    expect_error(improvedCV2(X[,0], spike.type="Spikes"), "need two or more cells to compute variances")

    out <- improvedCV2(counts, is.spike=!logical(ngenes))
    expect_true(all(is.na(out$p.value)))
    expect_error(improvedCV2(counts, is.spike=logical(ngenes)), "need at least 2 spike-ins for trend fitting")
})

# Testing for robustness to zeroes.
test_that("improvedCV2 is robust to zeroes", {
    counts[1,] <- 0
    counts[501,] <- 0
    out <- improvedCV2(counts, is.spike=NA)
    expect_identical(out$mean[1], 0)
    expect_identical(out$var[1], 0)
    expect_identical(out$cv2[1], NaN)
    expect_identical(out$trend[1], NA_real_)
    expect_equivalent(out[1,], out[501,])
})
