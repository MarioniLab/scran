#include "scran.h"
#include "shuffle_custom.h"

template<class V, class I, class O>
void shuffle_matrix_internal(I in, O out, SEXP seed) {
    const size_t NR=in->get_nrow(), NC=in->get_ncol();
    V tmp(NR);

    boost::random::mt19937 gen(check_integer_scalar(seed, "seed"));
    for (size_t c=0; c<NC; ++c) {
        in->get_col(c, tmp.begin());
        shuffle_custom(tmp.begin(), tmp.end(), gen);
        out->set_col(c, tmp.begin());
    }

    return;
}

SEXP shuffle_matrix(SEXP incoming, SEXP seed) {
    BEGIN_RCPP
    int rtype=beachmat::find_sexp_type(incoming);
    if (rtype==INTSXP) {
        auto mat=beachmat::create_integer_matrix(incoming);
        auto out=beachmat::create_integer_output(mat->get_nrow(), mat->get_ncol(), beachmat::output_param(mat->get_matrix_type(), true, true));
        shuffle_matrix_internal<Rcpp::IntegerVector>(mat.get(), out.get(), seed);
        return out->yield();
    } else {
        auto mat=beachmat::create_numeric_matrix(incoming);
        auto out=beachmat::create_numeric_output(mat->get_nrow(), mat->get_ncol(), beachmat::output_param(mat->get_matrix_type(), true, true));
        shuffle_matrix_internal<Rcpp::NumericVector>(mat.get(), out.get(), seed);
        return out->yield();
    }
    END_RCPP
}
