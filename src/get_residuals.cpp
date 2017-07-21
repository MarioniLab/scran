#include "scran.h"
#include "run_dormqr.h"

/* This function computes residuals in a nice and quick manner.
 * It takes a pre-computed QR matrix from qr(), and optionally
 * a subset of integer indices, and returns a matrix of residuals.
 */

SEXP get_residuals(SEXP exprs, SEXP qr, SEXP qraux, SEXP subset, SEXP lower_bound) {
    BEGIN_RCPP
    auto emat=beachmat::create_numeric_matrix(exprs);
    const size_t& ncells=emat->get_ncol();

    // Checking the subset vector.
    auto subout=check_subset_vector(subset, emat->get_nrow());
    const size_t slen=subout.size();
    
    // Checking the QR matrix.
    run_dormqr multQ1(qr, qraux, 'T');
    run_dormqr multQ2(qr, qraux, 'N');
    const int ncoefs=multQ1.get_ncoefs();

    // Checking the lower bound.
    Rcpp::NumericVector lb(lower_bound);
    if (lb.size()!=1) {
        throw std::runtime_error("lower bound should be a numeric scalar");
    }
    const double lbound=lb[0];
    const bool check_lower=R_FINITE(lbound);

    auto omat=beachmat::create_numeric_output(slen, ncells, beachmat::output_param(exprs, true, false));
    Rcpp::NumericVector tmp(ncells);
    double* tptr=(ncells ? &(tmp[0]) : NULL);
    std::deque<int> below_bound;

    auto sIt=subout.begin();
    for (size_t s=0; s<slen; ++s, ++sIt) {
        emat->get_row(*sIt, tmp.begin());
            
        // Identifying elements below the lower bound.
        if (check_lower) { 
            auto tIt=tmp.begin();
            for (int c=0; c<ncells; ++c, ++tIt) {
                if (*tIt <= lbound) {
                    below_bound.push_back(c);                        
                }
            }
        }

        multQ1.run(tptr); // Getting main+residual effects.
        std::fill(tmp.begin(), tmp.begin()+ncoefs, 0); // setting main effects to zero.
        multQ2.run(tptr); // Getting residuals.

        // Forcing the values below the boundary to a value below the smallest residual.
        if (check_lower && !below_bound.empty()) {
            const double lowest=*std::min_element(tmp.begin(), tmp.end()) - 1;
            for (auto bbIt=below_bound.begin(); bbIt!=below_bound.end(); ++bbIt) { 
                tmp[*bbIt]=lowest;
            }
            below_bound.clear();
        }

        omat->set_row(s, tmp.begin());
    }

    return omat->yield();
    END_RCPP
}

