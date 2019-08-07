# Mocking up data for use throughout the various examples.
# Note that we use '::' as Roxygen seems to refuse to import
# anything when it encounters this script.
example.sce <- (function() {
    # Mocking up the count data.
    set.seed(100)
    ncells <- 200

    nspikes <- 100
    spike.means <- 2^stats::runif(nspikes, 3, 8)
    spike.disp <- 100/spike.means + 0.5
    spike.data <- matrix(stats::rnbinom(nspikes*ncells, mu=spike.means, size=1/spike.disp), ncol=ncells)

    ngenes <- 2000
    cell.means <- 2^stats::runif(ngenes, 2, 10)
    cell.disp <- 100/cell.means + 0.5
    cell.data <- matrix(stats::rnbinom(ngenes*ncells, mu=cell.means, size=1/cell.disp), ncol=ncells)

    combined <- rbind(cell.data, spike.data)
    colnames(combined) <- seq_len(ncells)
    rownames(combined) <- seq_len(nrow(combined))

    # Creating the SingleCellExperiment (minor hack to avoid the warning,
    # given that it doesn't really matter).
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts=combined))
    sce <- SingleCellExperiment::splitSCEByAlt(sce, 
        ifelse(seq_len(nrow(sce)) > ngenes, "Spike", "genes"))
    scater::logNormCounts(sce)
})()
