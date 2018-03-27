#include "scran.h"

template<class V, class I, class O>
void shuffle_matrix_internal(I in, O out) {
    const size_t NR=in->get_nrow(), NC=in->get_ncol();
    V tmp(NR);
    Rcpp::RNGScope rng; // Place after initialization of all Rcpp objects.

    for (size_t c=0; c<NC; ++c) {
        in->get_col(c, tmp.begin());
        Rx_shuffle(tmp.begin(), tmp.end());
        out->set_col(c, tmp.begin());
    }

    return;
}

SEXP shuffle_matrix(SEXP incoming) {
    BEGIN_RCPP
    int rtype=beachmat::find_sexp_type(incoming);
    if (rtype==INTSXP) {
        auto mat=beachmat::create_integer_matrix(incoming);
        auto out=beachmat::create_integer_output(mat->get_nrow(), mat->get_ncol(), beachmat::output_param(mat->get_matrix_type(), true, true));
        shuffle_matrix_internal<Rcpp::IntegerVector>(mat.get(), out.get());
        return out->yield();
    } else {
        auto mat=beachmat::create_numeric_matrix(incoming);
        auto out=beachmat::create_numeric_output(mat->get_nrow(), mat->get_ncol(), beachmat::output_param(mat->get_matrix_type(), true, true));
        shuffle_matrix_internal<Rcpp::NumericVector>(mat.get(), out.get());
        return out->yield();
    }
    END_RCPP
}
