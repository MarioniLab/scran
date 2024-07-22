# This checks the denoisePCA function.
# require(scran); require(testthat); source("setup.R"); source("test-denoise-pca.R")

set.seed(70001)
test_that("denoisePCANumber works as expected", {
    v <- sort(runif(100))
    total <- sum(v)

    tech.var <- total * 0.8
    out <- denoisePCANumber(v, tech.var, total)
    expect_identical(out, length(v) - sum(cumsum(rev(v)) < tech.var))

    alt <- denoisePCANumber(head(v, out+5L), tech.var, total)
    expect_identical(out, alt)
    alt <- denoisePCANumber(head(v, out+1L), tech.var, total)
    expect_identical(out, alt)

    alt <- denoisePCANumber(head(v, out-1L), tech.var, total)
    expect_identical(out-1L, alt)

    tech.var <- 0
    out <- denoisePCANumber(v, tech.var, total)
    expect_identical(out, length(v))

    tech.var <- total
    out <- denoisePCANumber(v, tech.var, total)
    expect_identical(out, 1L)
}) 

################################
# Running tests for denoisePCA. This requires a mean-variance trend, 
# hence the somewhat complex set-up for the mock data.

set.seed(1000)
ngenes <- 1000
npops <- 10
ncells <- 200
means <- 2^runif(ngenes, -1, 10)
pops <- matrix(2^rnorm(npops * ngenes), ncol=npops) * means

in.pop <- sample(npops, ncells, replace=TRUE)
true.means <- pops[,in.pop,drop=FALSE]

dispersions <- 10/means + 0.2
counts <- matrix(rnbinom(ngenes*ncells, mu=true.means, size=1/dispersions), ncol=ncells)
rownames(counts) <- paste0("Gene", seq_len(ngenes))

nspikes <- 100
chosen <- sample(ngenes, nspikes)
spikes <- matrix(rnbinom(nspikes*ncells, mu=true.means[chosen], size=1/dispersions[chosen]), ncol=ncells)
rownames(spikes) <- paste0("SPIKE", seq_len(nspikes))

# We use a spike-in-based trend to avoid potential issues with bio==0 when
# fitting directly to genes. These lead to fragile tests due to numerical
# imprecision breaking equality upon certain operations.
dec <- modelGeneVarWithSpikes(counts, spikes=spikes, 
    size.factors=rep(1, ncells), spike.size.factors=rep(1, ncells))
lcounts <- log2(counts + 1)

##########################################
##########################################

test_that("getDenoisedPCs works as expected", {
    d.out <- getDenoisedPCs(lcounts, technical=dec, subset.row=NULL)
    expect_identical(nrow(d.out$components), ncol(lcounts))

    verify_npcs <- function(d.out, sdev, tech.total) {
        npcs <- ncol(d.out$components)
        var.exp <- sdev^2
        total.var <- sum(var.exp)
        expect_equal(npcs[1], denoisePCANumber(var.exp, tech.total, total.var))

        # Chosen number of PCs should be at the technical threshold.
        expect_true(sum(var.exp[(npcs+1):ncol(lcounts)]) < tech.total) 
        expect_true(sum(var.exp[(npcs):ncol(lcounts)]) > tech.total)
    
        reported <- d.out$percent.var
        exp.var <- sdev^2
        expect_equal(reported, exp.var[seq_along(reported)]/sum(exp.var) * 100)
    }

    keep <- dec$bio > 0
    pc.out <- prcomp(t(lcounts[keep,]))
    total.tech <- sum(dec$tech[keep])
    verify_npcs(d.out, pc.out$sdev, total.tech)

    npcs <- ncol(d.out$components)
    expect_equal(d.out$components, pc.out$x[,seq_len(npcs)])
    expect_equivalent(d.out$rotation, pc.out$rotation[,seq_len(npcs)])

    # Checking with different values for the technical noise, just in case.
    for (sub in c(0.05, 0.1, 0.2)) {
        tmp <- dec
        tmp$tech <- tmp$tech - sub
        d.out2 <- getDenoisedPCs(lcounts, technical=tmp, subset.row=NULL)
        expect_false(ncol(d.out$components)==ncol(d.out2$components))

        keep <- tmp$total > tmp$tech
        pc.out <- prcomp(t(lcounts[keep,]))
        total.tech <- sum(tmp$tech[keep])
        verify_npcs(d.out2, pc.out$sdev, total.tech)

        npcs <- ncol(d.out2$components)
        expect_equal(d.out2$components, pc.out$x[,seq_len(npcs)])
        expect_equivalent(d.out2$rotation, pc.out$rotation[,seq_len(npcs)])
    }
})

test_that("Rotation vectors are projected correctly", {
    lrout <- getDenoisedPCs(lcounts, technical=metadata(dec)$trend, subset.row=NULL, fill.missing=TRUE)

    lcounts.extra <- rbind(lcounts, lcounts[1:10,])
    lrout.extra <- getDenoisedPCs(lcounts.extra, technical=metadata(dec)$trend, 
        subset.row=seq_len(nrow(lcounts)), fill.missing=TRUE)

    expect_equal(lrout$rotation[,], lrout.extra$rotation[seq_len(nrow(lcounts)),])
    expect_equal(lrout$rotation[1:10,], lrout.extra$rotation[nrow(lcounts)+seq_len(10),])

    # Checking that we get the exact input back when we ask for everything.
    lrout <- getDenoisedPCs(lcounts, technical=dec, min.rank=ncol(lcounts), 
        max.rank=ncol(lcounts), fill.missing=TRUE, subset.row=NULL)
    expect_equal(tcrossprod(lrout$rotation, lrout$components), lcounts - rowMeans(lcounts))
}) 

set.seed(1001)
test_that("getDenoisedPCs works with different technical inputs", {
    ref <- getDenoisedPCs(lcounts, technical=dec, subset.row=NULL)
    pcs <- getDenoisedPCs(lcounts, technical=dec$tech, subset.row=NULL)
    expect_equal(ref, pcs)

    # Row sums in C++ have different precision from row sums in R on 32 bit,
    # as the latter uses long double and thus 80-bit precision. This seems
    # to be enough to change the mean, and thus the technical trend, and thus
    # whether or not a gene is kept or retained. Insane stuff.
    #
    # Maybe this was fixed by my use of modelGeneVar above, but I'm not sure.
    if (.Platform$r_arch=="") {
        alt <- getDenoisedPCs(lcounts, technical=metadata(dec)$trend, subset.row=NULL)
        expect_equal(ref, alt)
    }

    # Testing the rescaling to force the total variance in dec$total to match the observed variance.
    # This occasionally fails if you are unfortunate to get something where tech==bio,
    # and the equality is broken when you do the rescaling manually.
    #
    # Maybe this was fixed by my use of modelGeneVar above, but I'm not willing to take the chance,
    # what with us being so close to release. 
    rescaled <- runif(nrow(lcounts))
    lcountsX <- lcounts * rescaled
    ref <- scran:::.get_denoised_pcs(lcountsX, technical=dec$tech * rescaled^2, subset.row=NULL)
    pcs <- scran:::.get_denoised_pcs(lcountsX, technical=dec, subset.row=NULL)
    expect_equal(ref, pcs)

    # Handles all-zero rows with zero variance, where scaling would be undefined.
    lcountsAlt <- lcounts    
    lcountsAlt[1,] <- 0
    decAlt <- dec	
    decAlt$total[1] <- decAlt$tech[1] <- decAlt$bio[1] <- 0
    ref <- getDenoisedPCs(lcountsAlt, technical=decAlt$tech, subset.row=NULL)
    pcs <- getDenoisedPCs(lcountsAlt, technical=decAlt, subset.row=NULL)
    expect_equal(ref, pcs)

    # Handles cases where observed variance is zero but reported variance is not, e.g., after blocking.
    lcountsAlt[1,] <- runif(ncol(lcountsAlt))
    ref <- getDenoisedPCs(lcountsAlt[-1,], technical=decAlt$tech[-1], subset.row=NULL)
    pcs <- getDenoisedPCs(lcountsAlt, technical=decAlt, subset.row=NULL)
    expect_equal(ref$components, pcs$components)
})

test_that("getDenoisedPCs works with subsetting", {
    sub <- sample(ngenes, ngenes/2)
    pcs <- getDenoisedPCs(lcounts, technical=dec, subset.row=sub)
    pcs2 <- getDenoisedPCs(lcounts[sub,], technical=dec[sub,], subset.row=NULL)

    are_PCs_equal(pcs$components, pcs2$components)
    are_PCs_equal(pcs$rotation, pcs2$rotation)
    expect_equal(pcs$percent.var, pcs2$percent.var)

    # Works with different technical inputs.
    alt.pcs <- getDenoisedPCs(lcounts, technical=metadata(dec)$trend, subset.row=sub)
    expect_equal(alt.pcs, pcs)
    alt.pcs <- getDenoisedPCs(lcounts, technical=dec$tech, subset.row=sub)
    expect_equal(alt.pcs, pcs)
})

test_that("getDenoisedPCs works with min/max rank settings", {
    # Setting the min/max at around ncol(ref) to force it to a predictable number of pcs.
    ref <- getDenoisedPCs(lcounts, technical=dec, subset.row=NULL)$components
    pcs <- getDenoisedPCs(lcounts, technical=dec, min.rank=ncol(ref)+1, subset.row=NULL)$components
    expect_identical(ncol(pcs), ncol(ref)+1L)
    expect_identical(pcs[,seq_len(ncol(ref))], ref[,]) 

    pcs <- getDenoisedPCs(lcounts, technical=dec, max.rank=ncol(ref)-1, subset.row=NULL)$components
    expect_identical(ncol(pcs), ncol(ref)-1L)
    expect_identical(pcs[,], ref[,-ncol(ref)])

    # Stress-testing some gibberish min/max settings.
    pcs <- getDenoisedPCs(lcounts, technical=dec, min.rank=ncol(lcounts), max.rank=ncol(ref), subset.row=NULL)$components
    expect_identical(ncol(pcs), ncol(ref))
    pcs <- getDenoisedPCs(lcounts, technical=dec, min.rank=ncol(ref), max.rank=Inf, subset.row=NULL)$components
    expect_identical(ncol(pcs), ncol(ref))
    pcs <- getDenoisedPCs(lcounts, technical=dec, min.rank=0, max.rank=Inf, subset.row=NULL)$components
    expect_identical(ncol(pcs), ncol(ref))
})

test_that("denoisePCA throws errors correctly", {
    expect_error(getDenoisedPCs(lcounts[0,], dec, subset.row=NULL), "same rows")
    expect_error(getDenoisedPCs(lcounts[0,], dec$tech, subset.row=NULL), "same as")
    expect_error(getDenoisedPCs(lcounts[0,,drop=FALSE], dec[0,], subset.row=NULL), "a dimension is zero")
    expect_error(getDenoisedPCs(lcounts[,0], dec, subset.row=NULL), "no residual d.f. in any level")

    expect_warning(getDenoisedPCs(lcounts, dec), "subset.row")
})

##########################################
##########################################

test_that("denoisePCA works with SingleCellExperiment inputs", {
    X <- SingleCellExperiment(list(logcounts=lcounts))
    expect_warning(X2 <- denoisePCA(X, technical=dec), "subset.row")
    pcx <- reducedDim(X2, "PCA")
    rownames(pcx) <- NULL

    pcs <- getDenoisedPCs(lcounts, technical=dec, fill.missing=TRUE, subset.row=NULL)
    are_PCs_equal(pcx, pcs$components)
    expect_identical(attr(pcx, "percentVar"), pcs$percent.var)

    # Checking lowrank calculations.
    set.seed(10)
    X3 <- denoisePCA(X, technical=dec, value="lowrank", subset.row=NULL)
    pcx <- assay(X3, "lowrank")
    expect_equivalent(as.matrix(pcx), tcrossprod(pcs$rotation, pcs$components))

    set.seed(10)
    X3b <- denoisePCA(rbind(X, X[1:10,]), technical=rbind(dec, dec[1:10,]), subset.row=1:nrow(X), value="lowrank")
    pcxb <- assay(X3b, "lowrank")
    expect_equivalent(as.matrix(pcx), as.matrix(pcxb)[1:nrow(X),]) 
    expect_equivalent(as.matrix(pcx[1:10,]), as.matrix(pcxb)[nrow(X) + 1:10,]) 

    X4 <- denoisePCA(X, technical=dec, value="lowrank", subset.row=1:200)
    expect_identical(dim(X3), dim(X4))

    X5 <- denoisePCA(X, technical=dec, value="lowrank", subset.row=1:200, preserve.shape=FALSE)
    expect_identical(rownames(X5), rownames(X)[seq_len(nrow(X)) <= 200L & dec$bio > 0])
})
