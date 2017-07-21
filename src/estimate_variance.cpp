#include "scran.h"
#include "run_dormqr.h"

template<class M>
SEXP estimate_variance_internal (SEXP qr, SEXP qraux, M emat, SEXP subset) {
    // Setting up for Q-based multiplication.
    run_dormqr multQ(qr, qraux, 'T');
    const int ncoefs=multQ.get_ncoefs();
    const int ncells=multQ.get_nobs();

    if (ncells!=int(emat->get_ncol())) {
        throw std::runtime_error("number of rows of QR matrix not equal to number of cells");
    } else if (ncells==0) {
        throw std::runtime_error("cannot compute variance for zero cells");
    }
    auto subout=check_subset_vector(subset, emat->get_nrow());
    const size_t slen=subout.size();

    // Setting up output objects.
    Rcpp::NumericVector means(slen), vars(slen);
    auto mIt=means.begin(), vIt=vars.begin();
    Rcpp::NumericVector tmp(ncells);

    // Running through each gene and reporting its variance and mean.
    for (auto sIt=subout.begin(); sIt!=subout.end(); ++sIt, ++mIt, ++vIt) {
        emat->get_row(*sIt, tmp.begin());
        (*mIt)=std::accumulate(tmp.begin(), tmp.end(), 0.0)/ncells;

        multQ.run(&(tmp[0])); // Okay, due to zero check above.
        double& curvar=(*vIt);
        for (auto tIt=tmp.begin()+ncoefs; tIt!=tmp.end(); ++tIt) { // only using the residual effects.
            curvar += (*tIt) * (*tIt);
        }
        curvar /= ncells - ncoefs;
    }

    return Rcpp::List::create(means, vars);
}

SEXP estimate_variance (SEXP qr, SEXP qraux, SEXP exprs, SEXP subset) {
    BEGIN_RCPP
    int rtype=beachmat::find_sexp_type(exprs);
    if (rtype==INTSXP) {
        auto emat=beachmat::create_integer_matrix(exprs);
        return estimate_variance_internal(qr, qraux, emat.get(), subset);
    } else {
        auto emat=beachmat::create_numeric_matrix(exprs);
        return estimate_variance_internal(qr, qraux, emat.get(), subset);
    }
    END_RCPP
}
