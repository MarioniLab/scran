# This checks the overlapExprs function.
# require(scran); require(testthat); source("test-overlap.R")

# Constructing a reference computer.
refcomp <- function(x, g1, g2, tol=1e-8) { 
    ngenes <- nrow(x)
    output <- numeric(ngenes)
    for (i in seq_len(ngenes)) {
        e1 <- x[i,g1]
        e2 <- x[i,g2]
        de <- outer(e1, e2, `-`)
        neq <- sum(abs(de) <= tol)
        ngr <- sum(de > tol)
        output[i] <- neq*0.5 + ngr
    }
    output <- output / (length(g1) * length(g2))
    names(output) <- rownames(x)
    return(output)
}

set.seed(11000)
Ngenes <- 20
Ncells <- 100
X <- log(matrix(rpois(Ngenes*Ncells, lambda=10), nrow=Ngenes)+1)
rownames(X) <- paste0("X", seq_len(Ngenes))
grouping <- as.character(sample(3, Ncells, replace=TRUE))

out <- scran:::.overlapExprs(X, grouping)
for (i1 in names(out)) { 
    obs <- out[[i1]]
    expect_identical(sort(c(i1, colnames(obs))), sort(unique(grouping)))
    for (i2 in colnames(obs)) { 
        ref <- refcomp(X, which(grouping==i1), which(grouping==i2))       
        expect_equal(ref, obs[,i2])
        expect_true(all(abs(out[[i1]][,i2] + out[[i2]][,i1] - 1) < 1e-8))
    }
}

# Checking subsetting works as expected.
chosen <- 17:11
alt <- scran:::.overlapExprs(X, grouping, subset.row=chosen)
out2 <- lapply(out, function(x) { x[chosen,] })
expect_equal(alt, out2)

#############################
# Checking what happens with a blocking factor.

blockcomp <- function(X, groups, block) { 
    unique.groups <- sort(unique(groups))
    ngroups <- length(unique(groups))
    output <- used.cells <- vector("list", ngroups)
    ngenes <- nrow(X)
    for (g in seq_len(ngroups)) { 
        temp <- matrix(0, ngenes, ngroups-1L) 
        colnames(temp) <- unique.groups[-g]
        rownames(temp) <- rownames(X)
        output[[g]] <- temp
        temp.n <- integer(ngroups-1)
        names(temp.n) <- colnames(temp)
        used.cells[[g]] <- temp.n
    }
    names(output) <- names(used.cells) <- unique.groups

    for (b in levels(block)) { 
        current <- block==b
        cur.X <- X[,current,drop=FALSE]
        cur.groups <- groups[current]
        out <- scran:::.overlapExprs(cur.X, cur.groups) # Should be the same as running it in each block separately.
        
        for (g1 in names(out)) { 
            for (g2 in colnames(out[[g1]])) { 
                N <- sum(cur.groups==g1) + sum(cur.groups==g2)
                output[[g1]][,g2] <- output[[g1]][,g2] + out[[g1]][,g2] * N
                used.cells[[g1]][[g2]] <- used.cells[[g1]][[g2]] + N
            }
        }
    }

    for (g1 in names(output)) { 
        output[[g1]] <- t(t(output[[g1]])/used.cells[[g1]])
    }
    for (i1 in names(output)) { 
        for (i2 in colnames(output[[i1]])) {
            left <- output[[i1]][,i2]
            right <- output[[i2]][,i1]
            if (!any(is.na(left))) { 
                expect_true(all(abs(output[[i1]][,i2] + output[[i2]][,i1] - 1) < 1e-8))
            } else {
                expect_true(all(is.na(left)))
                expect_true(all(is.na(right)))
            }
        }
    }
    output
}

grouping <- rep(1:4, each=25)
block <- factor(rep(rep(LETTERS[1:4], c(2, 5, 8, 10)), 4))
X <- log(matrix(rpois(Ngenes*Ncells, lambda=10), nrow=Ngenes)+1)
rownames(X) <- paste0("X", seq_len(Ngenes))
out <- scran:::.overlapExprs(X, grouping, design=model.matrix(~block))
ref <- blockcomp(X, grouping, block) 
expect_equal(out, ref)

grouping <- rep(1:4, each=25)
block <- factor(rep(rep(LETTERS[1:4], c(2, 5, 8, 10)), 4))
block[grouping==4] <- "A" # Checking what happens when a group is only present in one blocking level.
X <- log(matrix(rpois(Ngenes*Ncells, lambda=10), nrow=Ngenes)+1)
rownames(X) <- paste0("X", seq_len(Ngenes))
out <- scran:::.overlapExprs(X, grouping, design=model.matrix(~block))
ref <- blockcomp(X, grouping, block) 
expect_equal(out, ref)

grouping <- rep(1:4, each=25)
block <- factor(rep(rep(LETTERS[1:4], c(2, 5, 8, 10)), 4))
grouping[block=="A"] <- 1 # Checking what happens when one blocking level only contains one group.
X <- log(matrix(rpois(Ngenes*Ncells, lambda=10), nrow=Ngenes)+1)
rownames(X) <- paste0("X", seq_len(Ngenes))
out <- scran:::.overlapExprs(X, grouping, design=model.matrix(~block))
ref <- blockcomp(X, grouping, block) 
expect_equal(out, ref)

grouping <- rep(1:4, each=25)
block <- factor(rep(rep(LETTERS[1:4], c(2, 5, 8, 10)), 4))
grouping[block=="A"] <- 1 # Checking what happens when blocks are confounded with group.
block[grouping==1] <- "A"
X <- log(matrix(rpois(Ngenes*Ncells, lambda=10), nrow=Ngenes)+1)
rownames(X) <- paste0("X", seq_len(Ngenes))
out <- scran:::.overlapExprs(X, grouping, design=model.matrix(~block))
ref <- blockcomp(X, grouping, block) 
expect_equal(out, ref)
expect_true(all(is.na(out[[1]])))

#############################
# Checking what happens when you force it to use residuals.

grouping <- rep(1:4, each=25)
block <- factor(rep(rep(LETTERS[1:4], c(2, 5, 8, 10)), 4))
design <- model.matrix(~block)
X <- log(matrix(rpois(Ngenes*Ncells, lambda=10), nrow=Ngenes)+1)
rownames(X) <- paste0("X", seq_len(Ngenes))

out <- scran:::.overlapExprs(X, grouping, design=design, residuals=TRUE)
fit <- lm.fit(y=t(X),x=design)
ref <- scran:::.overlapExprs(t(fit$residuals), grouping)
expect_equal(out, ref)

# A more natural example with a covariate.
block <- runif(Ncells)
design <- model.matrix(~block)
X <- log(matrix(rpois(Ngenes*Ncells, lambda=10), nrow=Ngenes)+1)
rownames(X) <- paste0("X", seq_len(Ngenes))

out <- scran:::.overlapExprs(X, grouping, design=design)
fit <- lm.fit(y=t(X),x=design)
ref <- scran:::.overlapExprs(t(fit$residuals), grouping)
expect_equal(out, ref)

#############################
# Checking for consistent behaviour with SCEsets.

Y <- matrix(rpois(Ngenes*Ncells, lambda=10), nrow=Ngenes)+1
rownames(X) <- paste0("X", seq_len(Ngenes))
X2 <- newSCESet(countData=Y, logExprsOffset=1, lowerDetectionLimit=0)
grouping <- rep(1:4, each=25)
expect_equal(overlapExprs(exprs(X2), grouping), overlapExprs(X2, grouping)) 
expect_equal(overlapExprs(counts(X2), grouping), overlapExprs(X2, grouping, assay="counts")) 

X2 <- calculateQCMetrics(X2, list(MySpike=rbinom(Ngenes, 1, 0.6)==0L))
setSpike(X2) <- "MySpike"
expect_equal(overlapExprs(exprs(X2), grouping, subset.row=!isSpike(X2)), overlapExprs(X2, grouping)) 
expect_equal(overlapExprs(exprs(X2), grouping), overlapExprs(X2, grouping, get.spikes=TRUE)) 

#############################
# Silly examples.

out <- scran:::.overlapExprs(X, integer(Ncells))
expect_identical(names(out), "0")
expect_identical(ncol(out[[1]]), 0L)

out <- scran:::.overlapExprs(X, grouping, subset.row=integer(0))
expect_identical(names(out), as.character(1:4))
expect_identical(unname(sapply(out, nrow)), integer(length(out))) 
out2 <- scran:::.overlapExprs(X[0,], grouping)
expect_identical(out, out2)

expect_identical(length(scran:::.overlapExprs(X[,0], grouping[0])), 0L)
expect_error(scran:::.overlapExprs(X[,0], grouping), "length of 'groups' not equal to number of cells", fixed=TRUE)
expect_error(scran:::.overlapExprs(X, grouping, design=cbind(rep(1, 10))), "'nrow(design)' not equal to number of cells", fixed=TRUE)


