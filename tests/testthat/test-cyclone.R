# This checks the cyclone implementation against a reference R-based implementation.
# require(scran); require(testthat); source("setup.R"); source("test-cyclone.R")

# On Windows:
# - C++ build failures for Win32 and Win64.
# - Weird test failures on Win32 with the optimized score calculation.
#
# On Mac:
# - Weird test failures for the optimized score calculation.
skip_on_os(c("windows", "mac"))

####################################################################################################

scramble_vector_src <- '
#include "Rcpp.h"
#include "boost/range/algorithm.hpp"
#include "pcg_random.hpp"
#include "convert_seed.h"
#include <algorithm>

// [[Rcpp::depends(BH, dqrng)]]
// [[Rcpp::export(rng=FALSE)]]
Rcpp::RObject scramble_vector(Rcpp::NumericVector invec, int niters, Rcpp::IntegerVector seed, int stream) {
    const size_t N=invec.size();
    Rcpp::NumericMatrix outmat(N, niters);
    const double* source=invec.begin();
    double* oIt=outmat.begin();

    auto generator=pcg32(dqrng::convert_seed<uint64_t>(seed), stream);
    for (int i=0; i<niters; ++i) {
        auto outcol=outmat.column(i);
        std::copy(source, source+N, outcol.begin());
        boost::range::random_shuffle(outcol, generator);
        source=oIt;
        oIt+=N;
    }

    return outmat;
}'

Rcpp::sourceCpp(code=scramble_vector_src)

####################################################################################################

classif.single <- function(cell, markers,Nmin.couples) { 
    left <- cell[markers[,1]]
    right <- cell[markers[,2]]

    tot <- sum(left!=right)
    if (tot < Nmin.couples) { 
        return(NA) 
    }  

    sum(left > right)/tot
}

random.success <- function(cell, markers, N, Nmin, Nmin.couples, seed, stream) {  
    test <- classif.single(cell,markers,Nmin.couples) 
    if (is.na(test)) { return(NA) } 

    cell.random <- scramble_vector(cell, N, seed, stream)
    success <- apply(cell.random, 2, classif.single, markers=markers, Nmin.couples=Nmin.couples)
    success <- success[!is.na(success)]
    if (length(success) < Nmin) { return(NA) }

    mean(success<test)
}

refFUN <- function(x, pairs) {
    x <- as.matrix(x)
    storage.mode(x) <- "double"
    gene.names <- rownames(x)

    chosen.x <- list()
    for (p in names(pairs)) {
        curp <- pairs[[p]]
        m1 <- match(curp$first, gene.names)
        m2 <- match(curp$second, gene.names)
        keep <- !is.na(m1) & !is.na(m2)
        m1 <- m1[keep]
        m2 <- m2[keep]
        
        all.present <- sort(unique(c(m1, m2)))
        chosen.x[[p]] <- x[all.present,,drop=FALSE]
        pairs[[p]] <- data.frame(first=match(m1, all.present),
                                 second=match(m2, all.present))
    }

    N <- 1000L
    Nmin <- 100L
    Nmin.couples <- 50L

    scores <- list()
    for (phase in c("G1", "S", "G2M")) {
        cur.x <- chosen.x[[phase]]
        cur.pairs <- pairs[[phase]]
        cur.scores <- numeric(ncol(x))

        rng.state <- scran:::.setup_pcg_state(ncol(x))
        rand.seeds <- rng.state$seeds[[1]]
        rand.streams <- rng.state$streams[[1]]

        for (i in seq_along(cur.scores)) {
            cur.scores[i] <- random.success(cell=cur.x[,i], markers=cur.pairs, N=N, Nmin=Nmin, Nmin.couples=Nmin.couples, 
                    seed=rand.seeds[[i]], stream=rand.streams[[i]])
        }
        scores[[phase]] <- cur.scores
    }

    scores <- do.call(data.frame, scores)
    scores.normalised <- scores/rowSums(scores)

    phases <- rep("S", ncol(x))
    phases[scores$G1 >= 0.5] <- "G1"
    phases[scores$G2M >= 0.5 & scores$G2M > scores$G1] <- "G2M"

    list(phases=phases, scores=scores, normalized.scores=scores.normalised)
}

####################################################################################################

# Spawning training data.

all.names <- paste0("X", seq_len(500))
Ngenes <- length(all.names)
all.pairs <- combn(Ngenes, 2)
re.pairs <- data.frame(first=all.names[all.pairs[1,]], second=all.names[all.pairs[2,]])

set.seed(100)
markers <- list(G1=re.pairs[sample(nrow(re.pairs), 100),],
                 S=re.pairs[sample(nrow(re.pairs), 200),],
               G2M=re.pairs[sample(nrow(re.pairs), 500),])

Ncells <- 10
test_that("cyclone works correctly on various datatypes", {
    # No ties.          
    set.seed(1000)
    X <- matrix(rnorm(Ngenes*Ncells), ncol=Ncells)
    rownames(X) <- all.names
    
    set.seed(100)
    reference <- refFUN(X, markers)
    set.seed(100)
    observed <- cyclone(X, markers)
    
    expect_identical(reference$phases, observed$phases)
    expect_equal(reference$scores, observed$scores)
    expect_equal(reference$normalized.scores, observed$normalized.scores)

    # Count data.
    set.seed(1001)
    X <- matrix(rpois(Ngenes*Ncells, lambda=10), ncol=Ncells)
    rownames(X) <- all.names
    
    set.seed(100)
    reference <- refFUN(X, markers)
    set.seed(100)
    observed <- cyclone(X, markers)
    
    expect_identical(reference$phases, observed$phases)
    expect_equal(reference$scores, observed$scores)
    expect_equal(reference$normalized.scores, observed$normalized.scores)

    # Low counts to induce more ties.
    set.seed(1002)
    X <- matrix(rpois(Ngenes*Ncells, lambda=1), ncol=Ncells)
    rownames(X) <- all.names
    
    set.seed(100)
    reference <- refFUN(X, markers)
    set.seed(100)
    observed <- cyclone(X, markers)
    
    expect_identical(reference$phases, observed$phases)
    expect_equal(reference$scores, observed$scores)
    expect_equal(reference$normalized.scores, observed$normalized.scores)
    
    # Changing the names of the marker sets.
    re.markers <- markers
    names(re.markers) <- paste0("X", names(markers))
    set.seed(100)
    re.out <- cyclone(X, re.markers)
    expect_identical(re.out$phase, character(0))
    expect_identical(colnames(re.out$scores), names(re.markers))
    expect_equal(unname(re.out$scores), unname(observed$scores))
})

set.seed(1003)
test_that("Cyclone gives the same results regardless of the number of cores", {
    X <- matrix(rpois(Ngenes*Ncells, lambda=100), ncol=Ncells)
    rownames(X) <- all.names

    set.seed(200)
    ref <- cyclone(X, markers)

    BPPARAM <- safeBPParam(3) # Before set.seed, as safeBPParam changes the seed.
    set.seed(200)
    alt <- cyclone(X, markers, BPPARAM=BPPARAM)
    expect_identical(ref, alt)
})

set.seed(1004)
test_that("Cyclone also works on SingleCellExperiment objects", {
    X <- matrix(rpois(Ngenes*Ncells, lambda=100), ncol=Ncells)
    rownames(X) <- all.names

    X2 <- SingleCellExperiment(list(counts=X))
    X2 <- scuttle::logNormCounts(X2)

    set.seed(100)
    reference <- refFUN(X, markers)
    
    set.seed(100)
    observed1 <- cyclone(X, markers)
    expect_equal(reference, observed1)
   
    # Doesn't matter whether you use the counts or logcounts. 
    set.seed(100)
    observed2 <- cyclone(X2, markers, assay.type="logcounts")
    expect_equal(reference, observed2)
})

set.seed(1005)
test_that("cyclone behaves correctly without cells or markers", {
    X <- matrix(rpois(Ngenes*Ncells, lambda=100), ncol=Ncells)
    rownames(X) <- all.names

    # Sensible behaviour with no cells.
    out <- cyclone(X[,0], markers)
    expect_identical(out$phases, character(0))
    expect_identical(nrow(out$scores), 0L)
    expect_identical(colnames(out$scores), c("G1", "S", "G2M"))
    expect_identical(nrow(out$normalized.scores), 0L)
    expect_identical(colnames(out$normalized.scores), c("G1", "S", "G2M"))
    
    # Sensible behaviour with no markers.
    no.markers <- list(G1=re.pairs[0,],
                        S=re.pairs[0,],
                      G2M=re.pairs[0,])
    out <- cyclone(X, no.markers)
    expect_true(all(is.na(out$phases)))
    expect_identical(colnames(out$scores), c("G1", "S", "G2M"))
    expect_identical(nrow(out$scores), ncol(X))
    expect_true(all(is.na(out$scores)))
    expect_identical(colnames(out$normalized.scores), c("G1", "S", "G2M"))
    expect_identical(nrow(out$normalized.scores), ncol(X))
    expect_true(all(is.na(out$normalized.scores)))
})

