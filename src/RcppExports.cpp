// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// build_snn_rank
Rcpp::List build_snn_rank(Rcpp::IntegerMatrix neighbors);
RcppExport SEXP _scran_build_snn_rank(SEXP neighborsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type neighbors(neighborsSEXP);
    rcpp_result_gen = Rcpp::wrap(build_snn_rank(neighbors));
    return rcpp_result_gen;
END_RCPP
}
// build_snn_number
Rcpp::List build_snn_number(Rcpp::IntegerMatrix neighbors);
RcppExport SEXP _scran_build_snn_number(SEXP neighborsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type neighbors(neighborsSEXP);
    rcpp_result_gen = Rcpp::wrap(build_snn_number(neighbors));
    return rcpp_result_gen;
END_RCPP
}
// calc_log_count_stats
Rcpp::List calc_log_count_stats(Rcpp::NumericVector Means, Rcpp::NumericVector Sizes, double tol, double disp, double pseudo);
RcppExport SEXP _scran_calc_log_count_stats(SEXP MeansSEXP, SEXP SizesSEXP, SEXP tolSEXP, SEXP dispSEXP, SEXP pseudoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type Means(MeansSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type Sizes(SizesSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< double >::type disp(dispSEXP);
    Rcpp::traits::input_parameter< double >::type pseudo(pseudoSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_log_count_stats(Means, Sizes, tol, disp, pseudo));
    return rcpp_result_gen;
END_RCPP
}
// calc_log_sqdiff
Rcpp::List calc_log_sqdiff(Rcpp::NumericVector Means, Rcpp::NumericVector Sizes, double tol, double disp, double pseudo, Rcpp::NumericVector Constants);
RcppExport SEXP _scran_calc_log_sqdiff(SEXP MeansSEXP, SEXP SizesSEXP, SEXP tolSEXP, SEXP dispSEXP, SEXP pseudoSEXP, SEXP ConstantsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type Means(MeansSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type Sizes(SizesSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< double >::type disp(dispSEXP);
    Rcpp::traits::input_parameter< double >::type pseudo(pseudoSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type Constants(ConstantsSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_log_sqdiff(Means, Sizes, tol, disp, pseudo, Constants));
    return rcpp_result_gen;
END_RCPP
}
// combine_rho
Rcpp::List combine_rho(int Ngenes, Rcpp::IntegerVector first, Rcpp::IntegerVector second, Rcpp::NumericVector Rho, Rcpp::NumericVector Pval, Rcpp::LogicalVector Limited, Rcpp::IntegerVector Order);
RcppExport SEXP _scran_combine_rho(SEXP NgenesSEXP, SEXP firstSEXP, SEXP secondSEXP, SEXP RhoSEXP, SEXP PvalSEXP, SEXP LimitedSEXP, SEXP OrderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< int >::type Ngenes(NgenesSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type first(firstSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type second(secondSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type Rho(RhoSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type Pval(PvalSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalVector >::type Limited(LimitedSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type Order(OrderSEXP);
    rcpp_result_gen = Rcpp::wrap(combine_rho(Ngenes, first, second, Rho, Pval, Limited, Order));
    return rcpp_result_gen;
END_RCPP
}
// combine_simes
Rcpp::NumericVector combine_simes(Rcpp::List Pvals, bool logp);
RcppExport SEXP _scran_combine_simes(SEXP PvalsSEXP, SEXP logpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type Pvals(PvalsSEXP);
    Rcpp::traits::input_parameter< bool >::type logp(logpSEXP);
    rcpp_result_gen = Rcpp::wrap(combine_simes(Pvals, logp));
    return rcpp_result_gen;
END_RCPP
}
// compute_CV2
Rcpp::List compute_CV2(SEXP exprs, Rcpp::IntegerVector subset_row, SEXP size_factors, SEXP log_prior);
RcppExport SEXP _scran_compute_CV2(SEXP exprsSEXP, SEXP subset_rowSEXP, SEXP size_factorsSEXP, SEXP log_priorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type exprs(exprsSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type subset_row(subset_rowSEXP);
    Rcpp::traits::input_parameter< SEXP >::type size_factors(size_factorsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type log_prior(log_priorSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_CV2(exprs, subset_row, size_factors, log_prior));
    return rcpp_result_gen;
END_RCPP
}
// compute_blocked_stats_lognorm
Rcpp::List compute_blocked_stats_lognorm(Rcpp::List bygroup, SEXP inmat, Rcpp::NumericVector sf, double pseudo);
RcppExport SEXP _scran_compute_blocked_stats_lognorm(SEXP bygroupSEXP, SEXP inmatSEXP, SEXP sfSEXP, SEXP pseudoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type bygroup(bygroupSEXP);
    Rcpp::traits::input_parameter< SEXP >::type inmat(inmatSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type sf(sfSEXP);
    Rcpp::traits::input_parameter< double >::type pseudo(pseudoSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_blocked_stats_lognorm(bygroup, inmat, sf, pseudo));
    return rcpp_result_gen;
END_RCPP
}
// compute_residual_stats_lognorm
Rcpp::List compute_residual_stats_lognorm(Rcpp::RObject qr, Rcpp::RObject qraux, SEXP inmat, Rcpp::NumericVector sf, double pseudo);
RcppExport SEXP _scran_compute_residual_stats_lognorm(SEXP qrSEXP, SEXP qrauxSEXP, SEXP inmatSEXP, SEXP sfSEXP, SEXP pseudoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject >::type qr(qrSEXP);
    Rcpp::traits::input_parameter< Rcpp::RObject >::type qraux(qrauxSEXP);
    Rcpp::traits::input_parameter< SEXP >::type inmat(inmatSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type sf(sfSEXP);
    Rcpp::traits::input_parameter< double >::type pseudo(pseudoSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_residual_stats_lognorm(qr, qraux, inmat, sf, pseudo));
    return rcpp_result_gen;
END_RCPP
}
// compute_blocked_stats_norm
Rcpp::List compute_blocked_stats_norm(Rcpp::List bygroup, SEXP inmat, Rcpp::NumericVector sf);
RcppExport SEXP _scran_compute_blocked_stats_norm(SEXP bygroupSEXP, SEXP inmatSEXP, SEXP sfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type bygroup(bygroupSEXP);
    Rcpp::traits::input_parameter< SEXP >::type inmat(inmatSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type sf(sfSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_blocked_stats_norm(bygroup, inmat, sf));
    return rcpp_result_gen;
END_RCPP
}
// compute_blocked_stats_none
Rcpp::List compute_blocked_stats_none(Rcpp::List bygroup, SEXP inmat);
RcppExport SEXP _scran_compute_blocked_stats_none(SEXP bygroupSEXP, SEXP inmatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type bygroup(bygroupSEXP);
    Rcpp::traits::input_parameter< SEXP >::type inmat(inmatSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_blocked_stats_none(bygroup, inmat));
    return rcpp_result_gen;
END_RCPP
}
// compute_residual_stats_none
Rcpp::List compute_residual_stats_none(Rcpp::RObject qr, Rcpp::RObject qraux, SEXP inmat);
RcppExport SEXP _scran_compute_residual_stats_none(SEXP qrSEXP, SEXP qrauxSEXP, SEXP inmatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject >::type qr(qrSEXP);
    Rcpp::traits::input_parameter< Rcpp::RObject >::type qraux(qrauxSEXP);
    Rcpp::traits::input_parameter< SEXP >::type inmat(inmatSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_residual_stats_none(qr, qraux, inmat));
    return rcpp_result_gen;
END_RCPP
}
// get_null_rho
Rcpp::NumericVector get_null_rho(int Ncells, int Niters, Rcpp::List Seeds, Rcpp::IntegerVector Streams);
RcppExport SEXP _scran_get_null_rho(SEXP NcellsSEXP, SEXP NitersSEXP, SEXP SeedsSEXP, SEXP StreamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< int >::type Ncells(NcellsSEXP);
    Rcpp::traits::input_parameter< int >::type Niters(NitersSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type Seeds(SeedsSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type Streams(StreamsSEXP);
    rcpp_result_gen = Rcpp::wrap(get_null_rho(Ncells, Niters, Seeds, Streams));
    return rcpp_result_gen;
END_RCPP
}
// get_null_rho_design
Rcpp::NumericVector get_null_rho_design(SEXP qr, SEXP qraux, int Niters, Rcpp::List Seeds, Rcpp::IntegerVector Streams);
RcppExport SEXP _scran_get_null_rho_design(SEXP qrSEXP, SEXP qrauxSEXP, SEXP NitersSEXP, SEXP SeedsSEXP, SEXP StreamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type qr(qrSEXP);
    Rcpp::traits::input_parameter< SEXP >::type qraux(qrauxSEXP);
    Rcpp::traits::input_parameter< int >::type Niters(NitersSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type Seeds(SeedsSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type Streams(StreamsSEXP);
    rcpp_result_gen = Rcpp::wrap(get_null_rho_design(qr, qraux, Niters, Seeds, Streams));
    return rcpp_result_gen;
END_RCPP
}
// test_rnorm
Rcpp::NumericVector test_rnorm(int N, SEXP seed, int stream);
RcppExport SEXP _scran_test_rnorm(SEXP NSEXP, SEXP seedSEXP, SEXP streamSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< SEXP >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< int >::type stream(streamSEXP);
    rcpp_result_gen = Rcpp::wrap(test_rnorm(N, seed, stream));
    return rcpp_result_gen;
END_RCPP
}
// compute_rho_pairs
Rcpp::NumericVector compute_rho_pairs(Rcpp::IntegerVector gene1, Rcpp::IntegerVector gene2, Rcpp::NumericMatrix ranks);
RcppExport SEXP _scran_compute_rho_pairs(SEXP gene1SEXP, SEXP gene2SEXP, SEXP ranksSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type gene1(gene1SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type gene2(gene2SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type ranks(ranksSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_rho_pairs(gene1, gene2, ranks));
    return rcpp_result_gen;
END_RCPP
}
// cyclone_scores
Rcpp::NumericVector cyclone_scores(Rcpp::IntegerVector mycells, SEXP exprs, SEXP marker1, SEXP marker2, SEXP indices, SEXP iter, SEXP miniter, SEXP minpair, SEXP seeds, SEXP streams);
RcppExport SEXP _scran_cyclone_scores(SEXP mycellsSEXP, SEXP exprsSEXP, SEXP marker1SEXP, SEXP marker2SEXP, SEXP indicesSEXP, SEXP iterSEXP, SEXP miniterSEXP, SEXP minpairSEXP, SEXP seedsSEXP, SEXP streamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type mycells(mycellsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type exprs(exprsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type marker1(marker1SEXP);
    Rcpp::traits::input_parameter< SEXP >::type marker2(marker2SEXP);
    Rcpp::traits::input_parameter< SEXP >::type indices(indicesSEXP);
    Rcpp::traits::input_parameter< SEXP >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< SEXP >::type miniter(miniterSEXP);
    Rcpp::traits::input_parameter< SEXP >::type minpair(minpairSEXP);
    Rcpp::traits::input_parameter< SEXP >::type seeds(seedsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type streams(streamsSEXP);
    rcpp_result_gen = Rcpp::wrap(cyclone_scores(mycells, exprs, marker1, marker2, indices, iter, miniter, minpair, seeds, streams));
    return rcpp_result_gen;
END_RCPP
}
// fit_linear_model
Rcpp::RObject fit_linear_model(Rcpp::RObject qr, SEXP qraux, SEXP exprs, SEXP subset, SEXP get_coefs);
RcppExport SEXP _scran_fit_linear_model(SEXP qrSEXP, SEXP qrauxSEXP, SEXP exprsSEXP, SEXP subsetSEXP, SEXP get_coefsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject >::type qr(qrSEXP);
    Rcpp::traits::input_parameter< SEXP >::type qraux(qrauxSEXP);
    Rcpp::traits::input_parameter< SEXP >::type exprs(exprsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type subset(subsetSEXP);
    Rcpp::traits::input_parameter< SEXP >::type get_coefs(get_coefsSEXP);
    rcpp_result_gen = Rcpp::wrap(fit_linear_model(qr, qraux, exprs, subset, get_coefs));
    return rcpp_result_gen;
END_RCPP
}
// fit_oneway
Rcpp::List fit_oneway(Rcpp::List grouping, SEXP exprs, SEXP subset);
RcppExport SEXP _scran_fit_oneway(SEXP groupingSEXP, SEXP exprsSEXP, SEXP subsetSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type grouping(groupingSEXP);
    Rcpp::traits::input_parameter< SEXP >::type exprs(exprsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type subset(subsetSEXP);
    rcpp_result_gen = Rcpp::wrap(fit_oneway(grouping, exprs, subset));
    return rcpp_result_gen;
END_RCPP
}
// get_residuals
Rcpp::RObject get_residuals(Rcpp::RObject exprs, SEXP qr, SEXP qraux, SEXP subset, SEXP lower_bound);
RcppExport SEXP _scran_get_residuals(SEXP exprsSEXP, SEXP qrSEXP, SEXP qrauxSEXP, SEXP subsetSEXP, SEXP lower_boundSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject >::type exprs(exprsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type qr(qrSEXP);
    Rcpp::traits::input_parameter< SEXP >::type qraux(qrauxSEXP);
    Rcpp::traits::input_parameter< SEXP >::type subset(subsetSEXP);
    Rcpp::traits::input_parameter< SEXP >::type lower_bound(lower_boundSEXP);
    rcpp_result_gen = Rcpp::wrap(get_residuals(exprs, qr, qraux, subset, lower_bound));
    return rcpp_result_gen;
END_RCPP
}
// get_scaled_ranks
Rcpp::RObject get_scaled_ranks(Rcpp::RObject exprs, Rcpp::RObject subset, Rcpp::RObject transpose, Rcpp::RObject as_sparse);
RcppExport SEXP _scran_get_scaled_ranks(SEXP exprsSEXP, SEXP subsetSEXP, SEXP transposeSEXP, SEXP as_sparseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject >::type exprs(exprsSEXP);
    Rcpp::traits::input_parameter< Rcpp::RObject >::type subset(subsetSEXP);
    Rcpp::traits::input_parameter< Rcpp::RObject >::type transpose(transposeSEXP);
    Rcpp::traits::input_parameter< Rcpp::RObject >::type as_sparse(as_sparseSEXP);
    rcpp_result_gen = Rcpp::wrap(get_scaled_ranks(exprs, subset, transpose, as_sparse));
    return rcpp_result_gen;
END_RCPP
}
// overlap_exprs
Rcpp::List overlap_exprs(Rcpp::RObject exprs, Rcpp::IntegerVector subset, Rcpp::List bygroup, Rcpp::RObject tolerance);
RcppExport SEXP _scran_overlap_exprs(SEXP exprsSEXP, SEXP subsetSEXP, SEXP bygroupSEXP, SEXP toleranceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject >::type exprs(exprsSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type subset(subsetSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type bygroup(bygroupSEXP);
    Rcpp::traits::input_parameter< Rcpp::RObject >::type tolerance(toleranceSEXP);
    rcpp_result_gen = Rcpp::wrap(overlap_exprs(exprs, subset, bygroup, tolerance));
    return rcpp_result_gen;
END_RCPP
}
// pool_size_factors
Rcpp::List pool_size_factors(Rcpp::RObject exprs, Rcpp::NumericVector pseudo_cell, Rcpp::IntegerVector order, Rcpp::IntegerVector pool_sizes);
RcppExport SEXP _scran_pool_size_factors(SEXP exprsSEXP, SEXP pseudo_cellSEXP, SEXP orderSEXP, SEXP pool_sizesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject >::type exprs(exprsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pseudo_cell(pseudo_cellSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type order(orderSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type pool_sizes(pool_sizesSEXP);
    rcpp_result_gen = Rcpp::wrap(pool_size_factors(exprs, pseudo_cell, order, pool_sizes));
    return rcpp_result_gen;
END_RCPP
}
// shuffle_matrix
Rcpp::RObject shuffle_matrix(Rcpp::RObject incoming, Rcpp::IntegerVector seed, int stream);
RcppExport SEXP _scran_shuffle_matrix(SEXP incomingSEXP, SEXP seedSEXP, SEXP streamSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject >::type incoming(incomingSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< int >::type stream(streamSEXP);
    rcpp_result_gen = Rcpp::wrap(shuffle_matrix(incoming, seed, stream));
    return rcpp_result_gen;
END_RCPP
}
// subset_and_divide
Rcpp::RObject subset_and_divide(Rcpp::RObject matrix, Rcpp::RObject row_subset, Rcpp::RObject col_subset, Rcpp::RObject scaling);
RcppExport SEXP _scran_subset_and_divide(SEXP matrixSEXP, SEXP row_subsetSEXP, SEXP col_subsetSEXP, SEXP scalingSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject >::type matrix(matrixSEXP);
    Rcpp::traits::input_parameter< Rcpp::RObject >::type row_subset(row_subsetSEXP);
    Rcpp::traits::input_parameter< Rcpp::RObject >::type col_subset(col_subsetSEXP);
    Rcpp::traits::input_parameter< Rcpp::RObject >::type scaling(scalingSEXP);
    rcpp_result_gen = Rcpp::wrap(subset_and_divide(matrix, row_subset, col_subset, scaling));
    return rcpp_result_gen;
END_RCPP
}
// test_shuffle_vector
Rcpp::RObject test_shuffle_vector(Rcpp::RObject incoming, Rcpp::RObject nits, Rcpp::RObject seed, Rcpp::RObject stream);
RcppExport SEXP _scran_test_shuffle_vector(SEXP incomingSEXP, SEXP nitsSEXP, SEXP seedSEXP, SEXP streamSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject >::type incoming(incomingSEXP);
    Rcpp::traits::input_parameter< Rcpp::RObject >::type nits(nitsSEXP);
    Rcpp::traits::input_parameter< Rcpp::RObject >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< Rcpp::RObject >::type stream(streamSEXP);
    rcpp_result_gen = Rcpp::wrap(test_shuffle_vector(incoming, nits, seed, stream));
    return rcpp_result_gen;
END_RCPP
}
// test_shuffle_matrix
Rcpp::RObject test_shuffle_matrix(Rcpp::RObject incoming, Rcpp::RObject seeds, Rcpp::RObject streams);
RcppExport SEXP _scran_test_shuffle_matrix(SEXP incomingSEXP, SEXP seedsSEXP, SEXP streamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject >::type incoming(incomingSEXP);
    Rcpp::traits::input_parameter< Rcpp::RObject >::type seeds(seedsSEXP);
    Rcpp::traits::input_parameter< Rcpp::RObject >::type streams(streamsSEXP);
    rcpp_result_gen = Rcpp::wrap(test_shuffle_matrix(incoming, seeds, streams));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_scran_build_snn_rank", (DL_FUNC) &_scran_build_snn_rank, 1},
    {"_scran_build_snn_number", (DL_FUNC) &_scran_build_snn_number, 1},
    {"_scran_calc_log_count_stats", (DL_FUNC) &_scran_calc_log_count_stats, 5},
    {"_scran_calc_log_sqdiff", (DL_FUNC) &_scran_calc_log_sqdiff, 6},
    {"_scran_combine_rho", (DL_FUNC) &_scran_combine_rho, 7},
    {"_scran_combine_simes", (DL_FUNC) &_scran_combine_simes, 2},
    {"_scran_compute_CV2", (DL_FUNC) &_scran_compute_CV2, 4},
    {"_scran_compute_blocked_stats_lognorm", (DL_FUNC) &_scran_compute_blocked_stats_lognorm, 4},
    {"_scran_compute_residual_stats_lognorm", (DL_FUNC) &_scran_compute_residual_stats_lognorm, 5},
    {"_scran_compute_blocked_stats_norm", (DL_FUNC) &_scran_compute_blocked_stats_norm, 3},
    {"_scran_compute_blocked_stats_none", (DL_FUNC) &_scran_compute_blocked_stats_none, 2},
    {"_scran_compute_residual_stats_none", (DL_FUNC) &_scran_compute_residual_stats_none, 3},
    {"_scran_get_null_rho", (DL_FUNC) &_scran_get_null_rho, 4},
    {"_scran_get_null_rho_design", (DL_FUNC) &_scran_get_null_rho_design, 5},
    {"_scran_test_rnorm", (DL_FUNC) &_scran_test_rnorm, 3},
    {"_scran_compute_rho_pairs", (DL_FUNC) &_scran_compute_rho_pairs, 3},
    {"_scran_cyclone_scores", (DL_FUNC) &_scran_cyclone_scores, 10},
    {"_scran_fit_linear_model", (DL_FUNC) &_scran_fit_linear_model, 5},
    {"_scran_fit_oneway", (DL_FUNC) &_scran_fit_oneway, 3},
    {"_scran_get_residuals", (DL_FUNC) &_scran_get_residuals, 5},
    {"_scran_get_scaled_ranks", (DL_FUNC) &_scran_get_scaled_ranks, 4},
    {"_scran_overlap_exprs", (DL_FUNC) &_scran_overlap_exprs, 4},
    {"_scran_pool_size_factors", (DL_FUNC) &_scran_pool_size_factors, 4},
    {"_scran_shuffle_matrix", (DL_FUNC) &_scran_shuffle_matrix, 3},
    {"_scran_subset_and_divide", (DL_FUNC) &_scran_subset_and_divide, 4},
    {"_scran_test_shuffle_vector", (DL_FUNC) &_scran_test_shuffle_vector, 4},
    {"_scran_test_shuffle_matrix", (DL_FUNC) &_scran_test_shuffle_matrix, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_scran(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}