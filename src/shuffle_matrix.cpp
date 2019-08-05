#include "Rcpp.h"

#include "beachmat/numeric_matrix.h"
#include "beachmat/integer_matrix.h"

#include "utils.h"
#include "rand_custom.h"

template<class I, class O>
Rcpp::RObject shuffle_matrix_internal(Rcpp::RObject incoming, Rcpp::IntegerVector seed, int stream) {
    auto in=beachmat::create_matrix<I>(incoming);
    const size_t NR=in->get_nrow(), NC=in->get_ncol();
    auto out=beachmat::create_output<O>(NR, NC, beachmat::output_param(in.get()));
    typename I::vector tmp(NR);

    auto gen=create_pcg32(seed, stream);
    for (size_t c=0; c<NC; ++c) {
        in->get_col(c, tmp.begin());
        shuffle_custom(tmp.begin(), tmp.end(), gen);
        out->set_col(c, tmp.begin());
    }

    return out->yield();
}

// [[Rcpp::export(rng=false)]]
Rcpp::RObject shuffle_matrix(Rcpp::RObject incoming, Rcpp::IntegerVector seed, int stream) {
    BEGIN_RCPP
    int rtype=beachmat::find_sexp_type(incoming);
    if (rtype==INTSXP) {
        return shuffle_matrix_internal<beachmat::integer_matrix, beachmat::integer_output>(incoming, seed, stream);
    } else {
        return shuffle_matrix_internal<beachmat::numeric_matrix, beachmat::numeric_output>(incoming, seed, stream);
    }
    END_RCPP
}
