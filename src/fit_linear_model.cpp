#include "Rcpp.h"

#include "beachmat/numeric_matrix.h"
#include "beachmat/integer_matrix.h"
#include "utils.h"
#include "run_dormqr.h"

#include <stdexcept>
#include <algorithm>
#include <vector>

template<class M>
SEXP fit_linear_model_internal (SEXP qr, SEXP qraux, SEXP inmat, SEXP subset, SEXP get_coefs) {
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
    auto subout=check_subset_vector(subset, emat->get_nrow());
    const size_t slen=subout.size();
    
    const bool coef_out=check_logical_scalar(get_coefs, "coefficient return specification");

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

// [[Rcpp::export(rng=false)]]
Rcpp::RObject fit_linear_model (Rcpp::RObject qr, SEXP qraux, SEXP exprs, SEXP subset, SEXP get_coefs) {
    int rtype=beachmat::find_sexp_type(exprs);
    if (rtype==INTSXP) {
        return fit_linear_model_internal<beachmat::integer_matrix>(qr, qraux, exprs, subset, get_coefs);
    } else {
        return fit_linear_model_internal<beachmat::numeric_matrix>(qr, qraux, exprs, subset, get_coefs);
    }
}

/* A much faster function when there's a one-way layout involved. */

template<class M>
Rcpp::List fit_oneway_internal (Rcpp::List bygroup, SEXP inmat, SEXP subset) {
    auto emat=beachmat::create_matrix<M>(inmat);
    const size_t ncells=emat->get_ncol();
 
    // Checking the various groupings.
    const size_t ngroups=bygroup.size();
    std::vector<Rcpp::IntegerVector> groups(ngroups);
    for (size_t i=0; i<ngroups; ++i) { 
        groups[i]=check_subset_vector(bygroup[i], ncells);
    }
    
    auto subout=check_subset_vector(subset, emat->get_nrow());
    const size_t slen=subout.size();
   
    // Setting up the output objects.
    Rcpp::NumericMatrix outvar(slen, ngroups);
    Rcpp::NumericMatrix outmean(slen, ngroups);
    int counter=0;
    Rcpp::NumericVector incoming(ncells);

    for (const auto& r : subout) {
        emat->get_row(r, incoming.begin());
        auto curvarrow=outvar.row(counter);
        auto curmeanrow=outmean.row(counter);
        ++counter;

        for (size_t g=0; g<ngroups; ++g) {
            const auto& curgroup=groups[g];
            double& curmean=curmeanrow[g];
            double& curvar=curvarrow[g];

            // Calculating the mean.          
            if (!curgroup.size()) {
                curmean=R_NaReal;
                curvar=R_NaReal;
                continue; 
            }
            for (const auto& index : curgroup) {
                curmean+=incoming[index];
            }
            curmean/=curgroup.size();

            // Computing the variance.
            if (curgroup.size()==1) {
                curvar=R_NaReal;
                continue;
            }
            for (const auto& index : curgroup) {
                const double tmp=incoming[index] - curmean;
                curvar += tmp * tmp;
            }
            curvar/=curgroup.size()-1;
        }
    }

    return Rcpp::List::create(outmean, outvar);
}

// [[Rcpp::export(rng=false)]]
Rcpp::List fit_oneway (Rcpp::List grouping, SEXP exprs, SEXP subset) {
    int rtype=beachmat::find_sexp_type(exprs);
    if (rtype==INTSXP) {
        return fit_oneway_internal<beachmat::integer_matrix>(grouping, exprs, subset);
    } else {
        return fit_oneway_internal<beachmat::numeric_matrix>(grouping, exprs, subset);
    }
}
