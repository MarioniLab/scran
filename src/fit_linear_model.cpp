#include "Rcpp.h"

#include "beachmat/numeric_matrix.h"
#include "beachmat/integer_matrix.h"
#include "utils.h"
#include "run_dormqr.h"

#include <stdexcept>
#include <algorithm>
#include <vector>

template<class M>
SEXP fit_linear_model_internal (SEXP qr, SEXP qraux, SEXP inmat, SEXP get_coefs) {
    // Setting up for Q-based multiplication.
    run_dormqr multQ(qr, qraux, 'T');
    const int ncoefs=multQ.get_ncoefs();
    const int ncells=multQ.get_nobs();

    auto emat=beachmat::create_matrix<M>(inmat);
    if (ncells!=static_cast<int>(emat->get_ncol())) {
        throw std::runtime_error("number of rows of QR matrix not equal to number of cells");
    } else if (ncells==0) {
        throw std::runtime_error("cannot compute variance for zero cells");
    }
    const size_t ngenes=emat->get_nrow();

    const bool coef_out=check_logical_scalar(get_coefs, "coefficient return specification");

    // Setting up output objects.
    Rcpp::NumericVector means(ngenes), vars(ngenes);
    auto mIt=means.begin(), vIt=vars.begin();
    Rcpp::NumericVector tmp(ncells);
    Rcpp::NumericMatrix coefs((coef_out ? ncoefs : 0), (coef_out ? ngenes : 0));
    auto cIt=coefs.begin();

    // Running through each gene and reporting its variance and mean.
    for (size_t s=0; s<ngenes; ++s) {
        emat->get_row(s, tmp.begin());
        (*mIt)=std::accumulate(tmp.begin(), tmp.end(), 0.0)/ncells;
        ++mIt;

        multQ.run(&(tmp[0])); // Okay, due to zero check above.
        double& curvar=(*vIt);
        for (auto tIt=tmp.begin()+ncoefs; tIt!=tmp.end(); ++tIt) { // only using the residual effects.
            curvar += (*tIt) * (*tIt);
        }
        curvar /= ncells - ncoefs;
        ++vIt;

        if (coef_out) { 
            multQ.solve(&(tmp[0]));
            std::copy(tmp.begin(), tmp.begin()+ncoefs, cIt);
            cIt += ncoefs;
        }
    }
    
    if (coef_out) {
        return Rcpp::List::create(coefs, means, vars);
    } else {
        return Rcpp::List::create(means, vars);
    }
}

// [[Rcpp::export(rng=false)]]
Rcpp::RObject fit_linear_model (Rcpp::RObject qr, SEXP qraux, SEXP exprs, SEXP get_coefs) {
    int rtype=beachmat::find_sexp_type(exprs);
    if (rtype==INTSXP) {
        return fit_linear_model_internal<beachmat::integer_matrix>(qr, qraux, exprs, get_coefs);
    } else {
        return fit_linear_model_internal<beachmat::numeric_matrix>(qr, qraux, exprs, get_coefs);
    }
}
