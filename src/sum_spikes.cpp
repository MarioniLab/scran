#include "scran.h"

template<class V, class M>
SEXP sum_spikes_internal(M mat, Rcpp::IntegerVector spikedex) {
    const size_t& ngenes = mat->get_nrow();
    const size_t& ncells = mat->get_ncol();
    int last=-1;
    for (const auto& curdex : spikedex) { 
        if (curdex <= last || curdex<0 || curdex>=ngenes) {
            throw std::runtime_error("'spikedex' should contain sorted indices in [0, ngenes)");
        } 
    }

    V incoming(ngenes);
    Rcpp::NumericVector output(ncells);
    auto oIt=output.begin();
    for (size_t c=0; c<ncells; ++c, ++oIt) {
        auto iIt=mat->get_const_col(c, incoming.begin());
        for (const auto& curdex : spikedex) {
            (*oIt)+=*(iIt + curdex);
        }
    }

    return output;
}

SEXP sum_spikes(SEXP counts, SEXP spikedex) {
    BEGIN_RCPP
    int rtype=beachmat::find_sexp_type(counts);
    if (rtype==INTSXP) {
        auto mat=beachmat::create_integer_matrix(counts);
        return sum_spikes_internal<Rcpp::IntegerVector>(mat.get(), spikedex);
    } else {
        auto mat=beachmat::create_numeric_matrix(counts);
        return sum_spikes_internal<Rcpp::NumericVector>(mat.get(), spikedex);
    }
    END_RCPP
}
