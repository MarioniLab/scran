# This tests the findMarkers function.

# require(scran); require(testthat); source("test-markers.R")

set.seed(70000)
ncells <- 200
ngenes <- 1000
means <- 2^runif(ngenes, -1, 5)
dummy <- matrix(rnbinom(ngenes*ncells, mu=means, size=5), ncol=ncells, nrow=ngenes)

rownames(dummy) <- paste0("X", seq_len(ngenes))
X <- newSCESet(countData=data.frame(dummy))
sizeFactors(X) <- colSums(dummy)
X <- normalize(X)

# Setting up a reference function.
library(limma)
REFFUN <- function(y, design, clust.vals, output, pval.type="any") { 
    lfit <- lmFit(y, design)
    for (host in clust.vals) {
        collected.lfc <- collected.p <- list()
        
        for (target in setdiff(clust.vals, host)) {
            con <- setNames(numeric(ncol(design)), colnames(design))
            con[[host]] <- 1
            con[[target]] <- -1
            
            fit2 <- contrasts.fit(lfit, con)
            fit2 <- eBayes(fit2, trend=TRUE, robust=TRUE)
            res <- topTable(fit2, n=Inf, sort.by="none")
            
            collected.lfc[[paste0("logFC.", target)]] <- res$logFC
            collected.p[[target]] <- res$P.Value
        }
        
        # Checking the value of 'top'.
        all.ranks <- lapply(collected.p, rank, ties.method="first")
        cur.out <- output[[host]]
        for (i in seq(1, ngenes, length.out=51)) {
            current <- rownames(X)[unique(unlist(lapply(all.ranks, FUN=function(x) { which(x <= i) })))]
            expect_identical(sort(current), sort(cur.out$Gene[cur.out$Top <= i]))
        }
        reordered <- order(do.call(pmin, all.ranks))
        
        # Checking the log-fold changes.
        for (target in names(collected.lfc)) {
            expect_equal(collected.lfc[[target]][reordered], cur.out[[target]])
        }
        
        # Checking the FDR, after Simes' or with the IUT.
        combined.p <- do.call(cbind, collected.p)
        if (pval.type=="any") { 
            pval <- apply(combined.p, 1, FUN=function(x) { min(p.adjust(x, method="BH")) })
        } else {
            pval <- apply(combined.p, 1, FUN=max)
        }
        adj.p <- p.adjust(pval, method="BH")
        expect_equal(adj.p[reordered], cur.out$FDR)
    }  
    return(TRUE)
}

clust <- kmeans(t(exprs(X)), centers=3)
out <- findMarkers(X, clusters=clust$cluster)

clusters <- as.factor(clust$cluster)
design <- model.matrix(~0 + clusters)
colnames(design) <- levels(clusters)
REFFUN(exprs(X), design, levels(clusters), out)

# Checking that the IUT calculation is correct.

out <- findMarkers(X, clusters=clust$cluster, pval.type="all")
REFFUN(exprs(X), design, levels(clusters), out, pval.type="all")

# Checking how it behaves with a design matrix.

clust <- kmeans(t(exprs(X)), centers=3)
block <- factor(sample(2, ncol(X), replace=TRUE))
out.des <- findMarkers(X, clusters=clust$cluster, design=model.matrix(~block))

clusters <- as.factor(clust$cluster)
design <- model.matrix(~0 + clusters + block)
colnames(design) <- c(levels(clusters), "block2")
REFFUN(exprs(X), design, levels(clusters), out.des)

# Checking that the results are the same regardless of how I provide inputs.

out <- findMarkers(X, clusters=clust$cluster)
out2 <- findMarkers(exprs(X), clusters=clust$cluster)
expect_identical(out, out2)

out <- findMarkers(X, clusters=clust$cluster, subset.row=100:1)
out2 <- findMarkers(X[100:1,], clusters=clust$cluster)
expect_identical(out, out2)

X <- calculateQCMetrics(X, feature_controls=list(ERCC=1:100))
setSpike(X) <- "ERCC"
out <- findMarkers(X, clusters=clust$cluster)
out2 <- findMarkers(exprs(X)[-(1:100),], clusters=clust$cluster)
expect_identical(out, out2)

# Checking consistency upon silly inputs.

expect_error(findMarkers(exprs(X)[0,], clusters=clust$cluster), "incorrect number of subscripts on matrix")
expect_error(findMarkers(exprs(X)[,0], clusters=integer(0)), "contrasts can be applied only to factors with 2 or more levels")

