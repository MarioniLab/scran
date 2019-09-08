# Tests combineMarkers().
# library(scran); library(testthat); source("test-combine-markers.R")

set.seed(1000)

# Setting up some reference data.
ngroups <- 10
groups <- cbind(rep(seq_len(ngroups), each=ngroups), rep(seq_len(ngroups), ngroups)) 
groups <- groups[groups[,1]!=groups[,2],]    

output <- list()
genes <- paste0("GENE_", 1:100)
for (x in seq_len(nrow(groups))) {
    output[[x]] <- data.frame(p.value=runif(length(genes)), logFC=rnorm(length(genes)), row.names=genes)
}

#######################################################

test_that("combineMarkers is the same as a reference impementation", {
    # With Simes' method.
    comb.any <- combineMarkers(output, groups)
    expect_identical(names(comb.any), as.character(seq_len(ngroups)))

    # With the IUT.
    comb.all <- combineMarkers(output, groups, pval.type="all")
    expect_identical(names(comb.any), as.character(seq_len(ngroups)))

    # Comparing to reference calculations.
    for (x in as.character(seq_len(ngroups))) {
        current.i <- as.character(groups[,1]) == x
        current.targets <- groups[current.i,2]
        current.stats <- output[current.i]

        # P-value reference calculations. 
        collected.p <- lapply(current.stats, "[[", i="p.value")
        pmat <- do.call(rbind, collected.p)
    
        obs.any <- comb.any[[x]][genes,]
        any.p <- apply(pmat, 2, FUN=function(p) { min(p.adjust(p, method="BH")) }) 
        expect_equal(any.p, obs.any$p.value)
        expect_equal(p.adjust(any.p, method="BH"), obs.any$FDR)

        obs.all <- comb.all[[x]][genes,]
        all.p <- apply(pmat, 2, FUN=max)
        expect_equal(all.p, obs.all$p.value)
        expect_equal(p.adjust(all.p, method="BH"), obs.all$FDR)

        # Rank calculations.
        expect_identical(rownames(comb.all[[x]]), genes[order(all.p)])

        min.rank <- do.call(pmin, lapply(collected.p, rank))
        min.p <- do.call(pmin, collected.p)
        any.o <- order(min.rank, min.p)
        expect_identical(rownames(comb.any[[x]]), genes[any.o])
        expect_equal(comb.any[[x]]$Top, min.rank[any.o])

        # Effect size extraction.
        for (other in seq_along(current.targets)) {
            cureffect <- paste0("logFC.", current.targets[other])
            expect_equal(obs.any[[cureffect]], obs.all[[cureffect]])
            expect_equal(current.stats[[other]]$logFC, obs.all[[cureffect]])
        }
    }

    # More Simes with some.
    comb.some <- combineMarkers(output, groups, pval.type="some")
    for (x in as.character(seq_len(ngroups))) {
        x.any <- comb.any[[x]]
        x.some <- comb.some[[x]]
        expect_identical(x.any[order(x.any$p.value),-1], x.some)
    }
})

test_that("combineMarkers responds to the pairing input", {
    # Handles different ordering of pairs.
    s <- sample(nrow(groups))
    comb.any <- combineMarkers(output, groups)
    alt.any <- combineMarkers(output[s], groups[s,])
    expect_identical(names(alt.any), as.character(unique(groups[s,1])))

    comb.all <- combineMarkers(output, groups, pval.type="all")
    alt.all <- combineMarkers(output[s], groups[s,], pval.type="all")
    expect_identical(names(alt.all), as.character(unique(groups[s,1])))

    first <- as.character(groups[,1])
    second <- as.character(groups[,2])
    for (x in names(comb.any)) {
        old.order <- second[first==x]
        cur.any <- comb.any[[x]]
        expect_equal(tail(colnames(cur.any), -3), paste0("logFC.", old.order))
        cur.all <- comb.all[[x]]
        expect_equal(tail(colnames(cur.all), -2), paste0("logFC.", old.order))

        new.order <- second[s][first[s]==x]
        cur.any.alt <- alt.any[[x]]
        expect_equal(tail(colnames(cur.any.alt), -3), paste0("logFC.", new.order))
        cur.all.alt <- alt.all[[x]]
        expect_equal(tail(colnames(cur.all.alt), -2), paste0("logFC.", new.order))

        expect_equal(cur.any, cur.any.alt[,colnames(cur.any)])
        expect_equal(cur.all, cur.all.alt[,colnames(cur.all)])
    }

    # Handles direct character input.
    groups2 <- LETTERS[groups]
    dim(groups2) <- dim(groups)
    re.comb <- combineMarkers(output, groups2)
    expect_identical(names(re.comb), LETTERS[as.integer(names(comb.any))])

    for (x in seq_along(re.comb)) {
        cur1 <- comb.any[[x]]
        cur2 <- re.comb[[x]]

        g <- as.integer(sub("logFC.", "", tail(colnames(cur1), -3)))
        g2 <- paste0("logFC.", LETTERS[g])
        expect_equal(tail(colnames(cur2), -3), g2)

        colnames(cur2) <- colnames(cur1)
        expect_equal(cur1, cur2)
    }
})

test_that("combineMarkers works with log-transformations", {
    login <- lapply(output, FUN=function(x) { x$p.value <- log(x$p.value); x })

    for (ptype in c("all", "any")) {
        # Log-transforming input.
        comb <- combineMarkers(output, groups, pval.type=ptype)
        comb.log <- combineMarkers(login, groups, pval.type=ptype, log.p.in=TRUE, log.p.out=FALSE)
        expect_equal(comb, comb.log)

        # Log-transforming output.
        comb.logout <- combineMarkers(output, groups, pval.type=ptype, log.p.out=TRUE)
        ref <- lapply(comb, FUN=function(x) {
            x$p.value <- log(x$p.value)
            x$FDR <- log(x$FDR)
            ii <- match(c("p.value", "FDR"), colnames(x))
            colnames(x)[ii] <- paste0("log.", colnames(x)[ii])
            x
        })
        ref <- SimpleList(ref)
        expect_equal(ref, comb.logout)

        # Log-transformed input and output.
        comb.log.both <- combineMarkers(login, groups, pval.type=ptype, log.p.in=TRUE)
        expect_equal(ref, comb.log.both)
    }
})

test_that("combineMarkers works with renamed fields", {
    renamed <- lapply(output, FUN=function(x) { 
        ii <- match(c("p.value", "logFC"), colnames(x))
        colnames(x)[ii] <- c("PValue", "LFC")
        x
    })
    comb <- combineMarkers(output, groups)
    recomb <- combineMarkers(renamed, groups, pval.field="PValue", effect.field="LFC")

    recomb.back <- lapply(comb, FUN=function(x) {
        colnames(x) <- sub("LFC", "logFC", colnames(x))
        x
    })
    recomb.back <- SimpleList(recomb.back)
    expect_identical(comb, recomb.back)

    # Handles the output field.
    re.recomb <- combineMarkers(renamed, groups, pval.field="PValue", effect.field="LFC", output.field="logFC")
    expect_identical(comb, re.recomb)

    # Checking that it runs without problems in pathological cases of name clashes.
    crazy.names <- c("p.value", 2:ngroups)
    extreme <- combineMarkers(output, array(crazy.names[groups], dim(groups)), output.field="log", log.p.out=TRUE)
    expect_identical(sum(colnames(extreme[["2"]])=="log.p.value"), 2L)
})

test_that("combineMarkers correctly returns the full stats", {
    stuff <- combineMarkers(output, groups, full.stats=TRUE)
    all.groups <- as.character(seq_len(ngroups)) 

    for (x in all.groups) {
        for (y in setdiff(all.groups, x)) {
            current <- stuff[[x]][,paste0("stats.", y)]
            current <- as.data.frame(current)
            correspondence <- which(groups[,1]==x & groups[,2]==y)            
            expect_identical(current, output[[correspondence]][rownames(current),])
        }
    }
})

test_that("combineMarkers works with silly inputs", {
    expect_error(combineMarkers(output[1], groups[0,]), "must be equal")
    expect_identical(combineMarkers(output[0], groups[0,]), setNames(SimpleList(), character(0)))

    empty <- combineMarkers(lapply(output, FUN=function(x){ x[0,] }), groups)
    expect_identical(names(empty), as.character(seq_len(ngroups)))
    expect_identical(unname(vapply(empty, nrow, 0L)), integer(ngroups))
    expect_equal(unname(vapply(empty, ncol, 0L)), rep(ngroups + 2L, ngroups))

    empty.full <- combineMarkers(lapply(output, FUN=function(x){ x[0,] }), groups, full.stats=TRUE)
    expect_identical(names(empty), names(empty.full))
    expect_identical(unname(vapply(empty.full, nrow, 0L)), integer(ngroups))
    expect_equal(unname(vapply(empty.full, ncol, 0L)), rep(ngroups + 2L, ngroups))

    # Ignores missing levels entirely.
    as.df <- data.frame(factor(groups[,1], levels=seq_len(ngroups+1)), factor(groups[,2]))
    ref <- combineMarkers(output, groups)
    lost <- combineMarkers(output, as.df)
    expect_identical(ref, lost)

    # Handles NA values in the group specifier.
    groups0 <- groups
    groups0[1,] <- NA
    groups0[2,1] <- NA
    groups0[2,2] <- NA
    ref <- combineMarkers(output, groups0)
    lost <- combineMarkers(output[-(1:2)], groups0[-(1:2),])
    expect_identical(ref, lost)

    # Handles NA values in the statistics.
    output0 <- output
    output0[[1]]$p.value <- NA_real_
    ref <- combineMarkers(output0, groups)
    ref[["1"]]$logFC.2 <- NULL
    lost <- combineMarkers(output0[-1], groups[-1,])
    expect_identical(ref, lost)
})
