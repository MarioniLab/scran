#include "scran.h"
#include "run_dormqr.h"

template<class M>
SEXP fit_linear_model_internal (SEXP qr, SEXP qraux, M emat, SEXP subset, SEXP get_coefs) {
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
    
    Rcpp::LogicalVector gcoef(get_coefs);
    if (gcoef.size()!=1) {
        throw std::runtime_error("'get_coefs' should be a logical scalar");
    }
    const bool coef_out=gcoef[0];

    // Setting up output objects.
    Rcpp::NumericVector means(slen), vars(slen);
    auto mIt=means.begin(), vIt=vars.begin();
    Rcpp::NumericVector tmp(ncells);
    Rcpp::NumericMatrix coefs((coef_out ? ncoefs : 0), (coef_out ? slen : 0));
    auto cIt=coefs.begin();

    // Running through each gene and reporting its variance and mean.
    for (const auto& s : subout) { 
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

SEXP fit_linear_model (SEXP qr, SEXP qraux, SEXP exprs, SEXP subset, SEXP get_coefs) {
    BEGIN_RCPP
    int rtype=beachmat::find_sexp_type(exprs);
    if (rtype==INTSXP) {
        auto emat=beachmat::create_integer_matrix(exprs);
        return fit_linear_model_internal(qr, qraux, emat.get(), subset, get_coefs);
    } else {
        auto emat=beachmat::create_numeric_matrix(exprs);
        return fit_linear_model_internal(qr, qraux, emat.get(), subset, get_coefs);
    }
    END_RCPP
}
