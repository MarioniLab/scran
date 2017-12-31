# This checks the overlapExprs function.
# require(scran); require(testthat); source("test-overlap.R")

set.seed(11000)
Ngenes <- 20
Ncells <- 100
X <- log(matrix(rpois(Ngenes*Ncells, lambda=10), nrow=Ngenes)+1)
rownames(X) <- paste0("X", seq_len(Ngenes))
grouping <- as.character(sample(3, Ncells, replace=TRUE))

#############################

refcomp <- function(x, g1, g2, tol=1e-8) 
# Constructing a reference calculator for each pair of groups.
{ 
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
    return(output)
}

RENAME <- function(i) { paste0("overlap.", i) }

checkSymmetry <- function(out) 
# Checks symmetry of the overlap values.
{
    for (i1 in names(out)) { 
        for (i2 in names(out)) { 
            if (i1==i2) { next }
            left <- out[[i1]][rownames(X),RENAME(i2)] 
            right <- out[[i2]][rownames(X),RENAME(i1)]
            expect_identical(is.na(left), is.na(right))
            expect_true(all(abs(left + right - 1) < 1e-8 | is.na(left)))
        }
    }
    return(invisible(NULL))
}

test_that("overlapExprs works in the simple case", {
    out <- overlapExprs(X, grouping)
    
    # Checking that all values are as expected.
    checkSymmetry(out)
    for (i1 in names(out)) { 
        obs <- out[[i1]]
        expect_identical(colnames(obs), c("Top", RENAME(setdiff(names(out), i1))))

        for (i2 in names(out)) { 
            if (i1==i2) { next }
            ref <- refcomp(X, which(grouping==i1), which(grouping==i2))       
            expect_equal(ref, obs[rownames(X),RENAME(i2)])
        }
    }

    # Checking subsetting works as expected.
    chosen <- 17:11
    alt <- overlapExprs(X, grouping, subset.row=chosen)
    for (i in names(alt)) {
        current <- alt[[i]]
        expect_identical(sort(rownames(current)), sort(rownames(X)[chosen]))
        expect_equal(current[,-1], out[[i]][rownames(current),-1])
    }

    chosen <- c(3,1,9,12,20)
    alt <- overlapExprs(X, grouping, subset.row=chosen)
    for (i in names(alt)) {
        current <- alt[[i]]
        expect_identical(sort(rownames(current)), sort(rownames(X)[chosen]))
        expect_equal(current[,-1], out[[i]][rownames(current),-1])
    }
})

#############################

checkValues <- function(out, ref, ignore.first=FALSE) 
# Checks identity of two sets of values, eliminating the difference in row orders.
{
    for (i1 in names(out)) { 
        obs <- out[[i1]]
        if (ignore.first) {
            expect_identical(colnames(obs)[-1], colnames(ref[[i1]])[-1])
        } else {
            expect_identical(colnames(obs), colnames(ref[[i1]]))
        }
        expect_equal(obs[,-1], ref[[i1]][rownames(obs),-1])        
    }
    return(invisible(NULL))
}

checkRanking <- function(out, FUN) 
# Checks the rankings for rank.type="any".
{ 
    for (i1 in names(out)) { 
        obs <- out[[i1]]
        top.rank <- best.val <- rep(NA_real_, nrow(obs))

        for (i2 in names(out)) {
            if (i1==i2) { next }
            curval <- FUN(obs[rownames(X),RENAME(i2)])
            top.rank <- pmin(top.rank, rank(curval, ties.method="first", na.last="keep"), na.rm=TRUE)
            best.val <- pmin(best.val, curval, na.rm=TRUE) 
        }

        o <- order(top.rank, best.val)
        expect_identical(rownames(X)[o], rownames(obs))
        expect_identical(obs$Top, as.integer(top.rank[o]))
    }
    return(invisible(NULL))
}

defaultRank <- function(x) { 0.5 - abs(x-0.5) }

test_that("overlapExprs responds to the request type", {
    ref <- overlapExprs(X, grouping)
    checkRanking(ref, defaultRank)
    
    # Checking the direction going down.
    out <- overlapExprs(X, grouping, direction="down")
    checkValues(out, ref)
    checkRanking(out, identity)

    # Checking the direction going up.
    out <- overlapExprs(X, grouping, direction="up")
    checkValues(out, ref)
    checkRanking(out, function(x) { 1 - x })
    
    # Checking the direction with rank.type="all".
    out <- overlapExprs(X, grouping, rank.type="all")
    checkValues(out, ref, ignore.first=TRUE)  
    for (i1 in names(out)) {
        current <- out[[i1]]
        gunk <- apply(current[,-1], 1, FUN=function(x) { x[which.max(defaultRank(x))] })
        expect_equal(names(gunk), rownames(current))
        expect_equal(unname(gunk), current$Worst)
    }

    out <- overlapExprs(X, grouping, rank.type="all", direction="up")
    checkValues(out, ref, ignore.first=TRUE)  
    for (i1 in names(out)) {
        current <- out[[i1]]
        gunk <- apply(current[,-1], 1, FUN=function(x) { x[which.max(1-x)] })
        expect_equal(names(gunk), rownames(current))
        expect_equal(unname(gunk), current$Worst)
    }
})

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
        out <- overlapExprs(cur.X, cur.groups) # Should be the same as running it in each block separately.
        
        for (g1 in names(out)) { 
            for (g2 in names(out)) { 
                if (g1==g2) { next }
                N <- sum(cur.groups==g1) * sum(cur.groups==g2)
                
                # using 'rownames' to force it back to the original order.
                output[[g1]][,g2] <- output[[g1]][,g2] + out[[g1]][rownames(X),RENAME(g2)] * N 
                used.cells[[g1]][[g2]] <- used.cells[[g1]][[g2]] + N
            }
        }
    }

    # Coercing to the same type of output format.
    for (g1 in names(output)) { 
        output[[g1]] <- t(t(output[[g1]])/used.cells[[g1]])
        colnames(output[[g1]]) <- RENAME(colnames(output[[g1]]))
        output[[g1]] <- DataFrame(Top=NA_integer_, output[[g1]], row.names=rownames(X))
    }
    output
}

test_that("overlapExprs works correctly with a blocking factor", {
    grouping <- rep(1:4, each=25)
    block <- factor(rep(rep(LETTERS[1:4], c(2, 5, 8, 10)), 4))
    X <- log(matrix(rpois(Ngenes*Ncells, lambda=10), nrow=Ngenes)+1)
    rownames(X) <- paste0("X", seq_len(Ngenes))

    out <- overlapExprs(X, grouping, block=block)
    ref <- blockcomp(X, grouping, block) 
    checkValues(out, ref)
    checkSymmetry(out)
    checkRanking(out, defaultRank) 
  
    # Checking that the rankings are still operational here. 
    out2 <- overlapExprs(X, grouping, block=block, direction="up")
    checkValues(out, out2)
    checkRanking(out2, function(x) { 1 - x })
 
    # Checking what happens when a group is only present in one blocking level.
    grouping <- rep(1:4, each=25)
    block <- factor(rep(rep(LETTERS[1:4], c(2, 5, 8, 10)), 4))
    block[grouping==4] <- "A" 
    X <- log(matrix(rpois(Ngenes*Ncells, lambda=10), nrow=Ngenes)+1)
    rownames(X) <- paste0("X", seq_len(Ngenes))

    out <- overlapExprs(X, grouping, block=block)
    ref <- blockcomp(X, grouping, block) 
    checkValues(out, ref)
    checkSymmetry(out)
    checkRanking(out, defaultRank) 
    
    # Checking what happens when one blocking level only contains one group.
    grouping <- rep(1:4, each=25)
    block <- factor(rep(rep(LETTERS[1:4], c(2, 5, 8, 10)), 4))
    grouping[block=="A"] <- 1 
    X <- log(matrix(rpois(Ngenes*Ncells, lambda=10), nrow=Ngenes)+1)
    rownames(X) <- paste0("X", seq_len(Ngenes))

    out <- overlapExprs(X, grouping, block=block)
    ref <- blockcomp(X, grouping, block) 
    checkValues(out, ref)
    checkSymmetry(out)
    checkRanking(out, defaultRank) 
    
    # Checking what happens when blocks are confounded with group.
    grouping <- rep(1:4, each=25)
    block <- factor(rep(rep(LETTERS[1:4], c(2, 5, 8, 10)), 4))
    grouping[block=="A"] <- 1 
    block[grouping==1] <- "A"
    X <- log(matrix(rpois(Ngenes*Ncells, lambda=10), nrow=Ngenes)+1)
    rownames(X) <- paste0("X", seq_len(Ngenes))
    
    out <- overlapExprs(X, grouping, block=block)
    ref <- blockcomp(X, grouping, block) 
    checkValues(out, ref)
    checkSymmetry(out)
    checkRanking(out, defaultRank) 
    expect_true(all(is.na(out[[1]][,-1])))
})

#############################
# Checking what happens when you force it to use residuals.

test_that("overlapExprs works correctly with residuals", {
    grouping <- rep(1:4, each=25)
    block <- factor(rep(rep(LETTERS[1:4], c(2, 5, 8, 10)), 4))
    design <- model.matrix(~block)
    X <- log(matrix(rpois(Ngenes*Ncells, lambda=10), nrow=Ngenes)+1)
    rownames(X) <- paste0("X", seq_len(Ngenes))
    
    # Unbounded.
    out <- overlapExprs(X, grouping, design=design, lower.bound=NA)
    fit <- lm.fit(y=t(X),x=design)
    ref <- overlapExprs(t(fit$residuals), grouping)
    expect_equal(out, ref)
    
    # Lower-bounded.
    X[] <- log(matrix(rpois(Ngenes*Ncells, lambda=1), nrow=Ngenes)+1)
    out <- overlapExprs(X, grouping, design=design, lower.bound=0) 
    fit <- lm.fit(y=t(X),x=design)
    resid <- t(fit$residuals)
    resid[X<=0] <- -100
    ref <- overlapExprs(resid, grouping)
    expect_equal(out, ref)
    
    ## A more natural example with a covariate.
    block <- runif(Ncells)
    design <- model.matrix(~block)
    X <- log(matrix(rpois(Ngenes*Ncells, lambda=10), nrow=Ngenes)+1)
    rownames(X) <- paste0("X", seq_len(Ngenes))
    
    out <- overlapExprs(X, grouping, design=design, lower.bound=NA)
    fit <- lm.fit(y=t(X),x=design)
    ref <- overlapExprs(t(fit$residuals), grouping)
    expect_equal(out, ref)
    
    X[] <- log(matrix(rpois(Ngenes*Ncells, lambda=1), nrow=Ngenes)+1)
    out <- overlapExprs(X, grouping, design=design, lower.bound=0) # With boundedness.
    fit <- lm.fit(y=t(X),x=design)
    resid <- t(fit$residuals)
    resid[X<=0] <- -100
    ref <- overlapExprs(resid, grouping)
    expect_equal(out, ref)
})

#############################
# Checking for consistent behaviour with SingleCellExperiment objects. 

test_that("overlapExprs behaves consistently with SingleCellExperiment objects", {
    Y <- matrix(rpois(Ngenes*Ncells, lambda=10), nrow=Ngenes)+1
    rownames(Y) <- paste0("X", seq_len(Ngenes))
    X2 <- SingleCellExperiment(list(counts=Y))
    sizeFactors(X2) <- colSums(Y)
    X2 <- normalize(X2)

    grouping <- rep(1:4, each=25)
    expect_equal(overlapExprs(exprs(X2), grouping), overlapExprs(X2, grouping)) 
    expect_equal(overlapExprs(counts(X2), grouping), overlapExprs(X2, grouping, assay.type="counts")) 

    isSpike(X2, "MySpike") <- rbinom(Ngenes, 1, 0.6)==0L
    expect_equal(overlapExprs(exprs(X2), grouping, subset.row=!isSpike(X2)), overlapExprs(X2, grouping)) 
    expect_equal(overlapExprs(exprs(X2), grouping), overlapExprs(X2, grouping, get.spikes=TRUE)) 
})

#############################
# Silly examples.

test_that("overlapExprs fails correctly on silly examples", {
    # Single-group example.
    out <- overlapExprs(X, integer(Ncells))
    expect_identical(names(out), "0")
    expect_identical(colnames(out[[1]]), "Top")
   
    #  No genes.
    grouping <- rep(1:4, each=25)
    out <- overlapExprs(X, grouping, subset.row=integer(0))
    expect_identical(names(out), as.character(1:4))
    expect_identical(unname(sapply(out, nrow)), integer(length(out))) 
    out2 <- overlapExprs(X[0,], grouping)
    expect_identical(out, out2)
    
    # No cells, or mismatched numbers of cells.
    expect_identical(length(overlapExprs(X[,0], grouping[0])), 0L)
    expect_error(overlapExprs(X[,0], grouping), "length of 'groups' not equal to number of cells", fixed=TRUE)
    expect_error(overlapExprs(X, grouping, design=cbind(rep(1, 10))), "'nrow(design)' not equal to number of cells", fixed=TRUE)
})

