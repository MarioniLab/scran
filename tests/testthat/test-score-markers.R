# This provides some tests for scoreMarkers.
# library(testthat); library(scran); source("test-score-markers.R")

test_that(".uniquify_DataFrame works as expected", {
    df <- DataFrame(group=rep(1:5, 3), block=rep(LETTERS[1:3], each=5))
    chosen <- sample(nrow(df), 1000, replace=TRUE)
    expanded <- df[chosen,]

    u.out <- scran:::.uniquify_DataFrame(expanded)
    expect_identical(expanded, u.out$unique[u.out$id,])
    expect_identical(sort(u.out$unique), sort(df))

    gblock <- scran:::.group_block_combinations(expanded$group, expanded$block)
    expect_identical(gblock$combinations, u.out$unique)
    expect_identical(gblock$id, u.out$id)

    # Works for a single column.
    df <- DataFrame(group=sample(letters, 100, replace=TRUE))
    u.out <- scran:::.uniquify_DataFrame(df)

    expect_identical(df[,1], u.out$unique[u.out$id,])
    expect_false(anyDuplicated(u.out$unique))

    gblock <- scran:::.group_block_combinations(df$group, NULL)
    expect_identical(gblock$combinations, u.out$unique)
    expect_identical(gblock$id, u.out$id)
})

test_that(".reindex_for_comparisons works as expected", {
    uniq.combos <- DataFrame(group=rep(1:5, 3), block=rep(LETTERS[1:3], each=5))

    desired <- DataFrame(left=1, right=2)
    df <- scran:::.reindex_comparisons_for_combinations(uniq.combos, desired)
    expect_true(all(uniq.combos$group[df$left]==2L))
    expect_true(all(uniq.combos$group[df$right]==1L))
    expect_identical(uniq.combos$block[df$left], LETTERS[1:3])
    expect_identical(uniq.combos$block[df$right], LETTERS[1:3])

    # Order doesn't matter.
    desired <- DataFrame(left=2, right=1)
    df2 <- scran:::.reindex_comparisons_for_combinations(uniq.combos, desired)
    expect_identical(df, df2)

    # Deduplicates correctly.
    desired <- DataFrame(left=1:2, right=2:1)
    df2 <- scran:::.reindex_comparisons_for_combinations(uniq.combos, desired)
    expect_identical(df, df2)

    # A more complex example.
    desired <- DataFrame(left=c(1, 5, 3), right=c(4, 4, 4))
    df <- scran:::.reindex_comparisons_for_combinations(uniq.combos, desired)
    group.left <- uniq.combos$group[df$left]
    group.right <- uniq.combos$group[df$right]
    expect_true(all((group.left %in% desired$right & group.right %in% desired$left) | (group.left %in% desired$left & group.right %in% desired$right)))
    expect_identical(uniq.combos$block[df$left], uniq.combos$block[df$right])
    expect_identical(uniq.combos$block[df$left], rep(LETTERS[1:3], each=3))

    # Works when not all blocks have the comparison of interest.
    uniq2 <- uniq.combos[-1,]
    desired <- DataFrame(left=c(1, 5, 3), right=c(2, 4, 4))
    ref <- scran:::.reindex_comparisons_for_combinations(uniq.combos, desired)
    df <- scran:::.reindex_comparisons_for_combinations(uniq2, desired)
    expect_identical(ref[-1,1], df[,1]+1L)
    expect_identical(ref[-1,2], df[,2]+1L)

    # Or when none of the blocks have the comparison of interest.
    desired <- DataFrame(left=6, right=1)
    df <- scran:::.reindex_comparisons_for_combinations(uniq.combos, desired)
    expect_identical(ncol(df), 2L)
    expect_identical(nrow(df), 0L)

    # Works without 'block='.
    desired <- DataFrame(left=2, right=1)
    df <- scran:::.reindex_comparisons_for_combinations(uniq.combos[,1,drop=FALSE], desired)
    expect_identical(df$left, 2L)
    expect_identical(df$right, 1L)
        
    # Same result with character names.
    desired <- DataFrame(left=c(5, 3, 2), right=c(2, 1, 4))
    ref <- scran:::.reindex_comparisons_for_combinations(uniq.combos, desired)

    scramble <- sample(letters)
    uniq.char <- DataFrame(group=scramble[uniq.combos$group], block=match(uniq.combos$block, LETTERS))
    desired.char <- DataFrame(left=scramble[desired$left], right=scramble[desired$right])
    df <- scran:::.reindex_comparisons_for_combinations(uniq.char, desired.char)

    expect_identical(ref, df)
})

test_that(".mapply_bind works as expected", {
    A1 <- matrix(runif(runif(20)), ncol=5)
    A2 <- matrix(runif(runif(20)), ncol=5)
    B1 <- matrix(runif(runif(10)), ncol=5)
    B2 <- matrix(runif(runif(10)), ncol=5)
    C1 <- matrix(runif(runif(15)), ncol=5)
    C2 <- matrix(runif(runif(15)), ncol=5)

    output <- scran:::.mapply_bind(df.list=list(list(X=A1, Y=A2), list(B1, B2), list(C1, C2)), FUN=rbind)
    expect_identical(output, list(X=rbind(A1, B1, C1), Y=rbind(A2, B2, C2)))

    A1 <- matrix(runif(runif(20)), nrow=5)
    A2 <- matrix(runif(runif(20)), nrow=5)
    B1 <- matrix(runif(runif(10)), nrow=5)
    B2 <- matrix(runif(runif(10)), nrow=5)
    C1 <- matrix(runif(runif(15)), nrow=5)
    C2 <- matrix(runif(runif(15)), nrow=5)

    output <- scran:::.mapply_bind(df.list=list(list(X=A1, Y=A2), list(B1, B2), list(C1, C2)), FUN=cbind)
    expect_identical(output, list(X=cbind(A1, B1, C1), Y=cbind(A2, B2, C2)))
})

.to_list <- function(x) unname(as.list(x))

test_that(".identify_effects_to_average works correctly", {
    scramble <- paste0("GROUP_", sample(letters)[1:5])
    uniq.combos <- DataFrame(group=rep(scramble, 3), block=rep(LETTERS[1:3], each=5))

    desired <- DataFrame(left=scramble[-2], right=scramble[2])
    reindexed <- scran:::.reindex_comparisons_for_combinations(uniq.combos, desired)
    averaged <- scran:::.identify_effects_to_average(uniq.combos, reindexed)

    expect_identical(lengths(averaged$indices.to.average), setNames(rep(3L, 4L), 1:4))
    expect_identical(averaged$averaged.comparisons[-1,], desired[-1,])
    expect_identical(.to_list(averaged$averaged.comparisons[1,]), rev(.to_list(desired[1,])))
    
    # Make a batch-specific comparison.
    uniq.combos2 <- rbind(uniq.combos, DataFrame(group="FOO", block="A"))
    desired2 <- DataFrame(left=c("FOO", scramble[1]), right=scramble[2])
    reindexed2 <- scran:::.reindex_comparisons_for_combinations(uniq.combos2, desired2)
    averaged2 <- scran:::.identify_effects_to_average(uniq.combos2, reindexed2)

    expect_identical(averaged2$averaged.comparisons[1,], desired2[1,])
    expect_identical(.to_list(averaged2$averaged.comparisons[-1,]), rev(.to_list(desired2[-1,])))
    expect_identical(averaged2$indices.to.average[[1]], 1L)
    expect_identical(averaged2$indices.to.average[[2]], 2:4)
})

test_that(".average_effects_across_blocks works correctly", {
    effects <- matrix(rnorm(1000), ncol=10)
    weights <- sample(20, 10)
    indices <- list(c(1, 5), c(6, 2, 7), c(4, 8))

    averaged <- scran:::.average_effect_across_blocks(effects, weights, indices)
    for (i in seq_along(indices)) {
        chosen <- indices[[i]]
        expect_equal(averaged$averaged.effects[[i]], DelayedMatrixStats::rowWeightedMeans(effects[,chosen], weights[chosen]))
        expect_equal(averaged$combined.weights[i], sum(weights[chosen]))
    }

    # Handles NA's in there correctly.
    effects[,1] <- NA
    effects[1:5,2] <- NA

    averaged <- scran:::.average_effect_across_blocks(effects, weights, indices)
    for (i in seq_along(indices)) {
        chosen <- indices[[i]]
        expect_equal(averaged$averaged.effects[[i]], DelayedMatrixStats::rowWeightedMeans(effects[,chosen], weights[chosen], na.rm=TRUE))
        expect_equal(averaged$combined.weights[i], sum(weights[chosen]))
    }
})

test_that(".cross_reference_to_desired works correctly", {
    desired <- DataFrame(left=c("A", "B", "C", "B"), right=c("D", "D", "D", "A"))
    observed <- DataFrame(left=c("A", "A", "D"), right=c("D", "B", "B"))

    out <- scran:::.cross_reference_to_desired(observed, desired)
    expect_identical(out$A$right, "D")
    expect_identical(out$A$direct.match, 1L)
    expect_identical(out$B$right, c("D", "A"))
    expect_identical(out$B$direct.match, c(NA_integer_, NA_integer_))
    expect_identical(out$B$flipped.match, 3:2)
    expect_identical(out$C$right, "D")
    expect_identical(out$C$direct.match, NA_integer_)
    expect_identical(out$C$flipped.match, NA_integer_)

    desired <- DataFrame(left=c("A", "B", "C", "B"), right=c("D", "A", "D", "D"))
    observed <- DataFrame(left=c("A", "B", "D"), right=c("D", "A", "B"))
    out <- scran:::.cross_reference_to_desired(observed, desired)
    expect_identical(out$B$right, c("A", "D"))
    expect_identical(out$B$direct.match, c(2L, NA_integer_))
    expect_identical(out$B$flipped.match, c(NA_integer_, 3L))
})

test_that(".cohen's D calculation works correctly", {
    means <- matrix(rnorm(100), ncol=5)
    vars <- matrix(runif(100), ncol=5)
    ncells <- sample(100, 5)

    left <- c(1,2,3,4,5)
    right <- c(2,3,4,5,1)
    left.ncells <- ncells[left]
    right.ncells <- ncells[right]

    out <- scran:::.compute_pairwise_cohen_d(means, vars, left, right, left.ncells, right.ncells)
    delta <- means[,left] - means[,right]
    pooled.s2 <- (t(vars[,left]) * (ncells[left] - 1) + t(vars[,right]) * (ncells[right] - 1)) /(ncells[left] + ncells[right] - 2)
    expect_identical(out, delta / t(sqrt(pooled.s2)))

    # Handles zero variances gracefully.
    means[1,1:2] <- 5
    vars[1,1:2] <- 0
    out <- scran:::.compute_pairwise_cohen_d(means, vars, left, right, left.ncells, right.ncells)
    expect_identical(out[1,1], 0)
    expect_false(any(out[-1]==0))

    means[1,1] <- 1
    out <- scran:::.compute_pairwise_cohen_d(means, vars, left, right, left.ncells, right.ncells)
    expect_identical(out[1,1], NaN)
    expect_false(any(is.na(out[-1])))
})

test_that("AUC calculations work correctly", {
    x <- matrix(rpois(1000, 0.5), ncol=100)
    combo.id <- sample(5, ncol(x), replace=TRUE)

    left <- c(5,4,3,2,1)
    right <- c(3,1,5,1,3)
    involved <- scran:::.group_by_used_combinations(combo.id, left, right, 5)

    auc <- scran:::.compute_auc(x, involved, left, right)
    for (i in seq_along(left)) {
        L <- involved[[left[i]]] + 1L
        R <- involved[[right[i]]] + 1L
        collected <- vapply(seq_len(nrow(x)), function(i) wilcox.test(x[i,L], x[i,R], exact=FALSE)$statistic, 0)
        expect_equal(collected / length(L) / length(R), auc[,i])
    }

    # Handles negative values.
    auc.neg <- scran:::.compute_auc(-x, involved, left, right)
    expect_equal(auc, 1 - auc.neg)

    # Handles not all groups being involved.
    chosen <- c(1, 3)
    involved <- scran:::.group_by_used_combinations(combo.id, left[chosen], right[chosen], 5)
    auc2 <- scran:::.compute_auc(x, involved, left[chosen], right[chosen])
    expect_identical(auc[,chosen], auc2)
})

test_that("logFC-detected calculations work correctly", {
    ncells <- sample(200, 5)
    ndetected <- t(matrix(rbinom(1000, ncells, 0.2), nrow=5))

    left <- c(5,4,3,2,1)
    right <- c(3,1,5,1,3)
    out <- scran:::.compute_lfc_detected(ndetected, left, right, ncells[left], ncells[right])

    ref.left <- t( t(ndetected[,left]) / ncells[left] )
    ref.right <- t( t(ndetected[,right]) / ncells[right] )

    # Shrinkage works correctly.
    pos <- out > 1e-8
    expect_identical(pos, ref.left > ref.right + 1e-8)
    expect_true(all(out[pos] < log2(ref.left[pos] / ref.right[pos])))

    neg <- out < -1e-8
    expect_identical(neg, ref.left < ref.right - 1e-8)
    expect_true(all(out[neg] > log2(ref.left[neg] / ref.right[neg])))

    # Properly set to zero.
    left <- c(5,3,1)
    right <- c(5,3,1)
    out <- scran:::.compute_lfc_detected(ndetected, left, right, ncells[left], ncells[right])
    expect_true(all(out==0))

    # Avoids undefined log-fold changes.
    copy <- ndetected
    copy[,1] <- 0
    left <- c(5,4,3,2,1)
    right <- c(3,1,5,1,3)
    out <- scran:::.compute_lfc_detected(copy, left, right, ncells[left], ncells[right])
    expect_false(any(!is.finite(out)))
})

test_that("overall effect size calculations work correctly", {
    y <- matrix(rpois(1000, 0.5), ncol=100)
    combo.id <- sample(5, ncol(x), replace=TRUE)
    ncells <- tabulate(combo.id, 5)

    left <- c(3,1,2)
    right <- c(5,5,1)
    left.ncells <- ncells[left]
    right.ncells <- ncells[right]

    stats <- scran:::compute_blocked_stats_none(y, combo.id - 1L, 5)
    cohen <- scran:::.compute_pairwise_cohen_d(stats[[1]], stats[[2]], left, right, left.ncells, right.ncells)
})
