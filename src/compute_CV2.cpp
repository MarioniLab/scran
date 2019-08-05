#include "Rcpp.h"

#include "beachmat/integer_matrix.h"
#include "beachmat/numeric_matrix.h"
#include "utils.h"

#include <stdexcept>
#include <cmath>
#include <algorithm>

/* This computes the mean and CV2 for every gene, while also dividing through by the size factors. 
 * Unless log_prior is specified, in which case the values in 'ptr' are assumed to be log-expressed.
 * These are converted back to the raw scale without any size factor calculations.
 */
template <class M>
Rcpp::List compute_CV2_internal(SEXP incoming, Rcpp::IntegerVector subset_row, SEXP size_factors, SEXP log_prior) {
    auto mat=beachmat::create_matrix<M>(incoming);
    auto rsubout=check_subset_vector(subset_row, mat->get_nrow());
    const size_t rslen=rsubout.size();
    const size_t ncells=mat->get_ncol();
    if (ncells < 2) {
        throw std::runtime_error("need two or more cells to compute variances");
    }

    // Checking the mode with which we will process the data.    
    const bool to_unlog=(log_prior!=R_NilValue);
    Rcpp::NumericVector sizefacs;
    double lp=0;
    if (to_unlog) { 
        lp=check_numeric_scalar(log_prior, "prior count");
        if (size_factors!=R_NilValue){ 
            throw std::runtime_error("size factors cannot be specified for log-expression input");
        }
    } else {
        sizefacs=size_factors;
        if (static_cast<size_t>(sizefacs.size())!=ncells) { 
            throw std::runtime_error("number of size factors is not equal to number of cells");
        }
    } 

    Rcpp::NumericVector means(rslen), vars(rslen);
    Rcpp::NumericVector tmp(ncells);
    auto mIt=means.begin(), vIt=vars.begin();
    
    for (auto rsIt=rsubout.begin(); rsIt!=rsubout.end(); ++rsIt, ++mIt, ++vIt) {
        mat->get_row(*rsIt, tmp.begin());

        if (!to_unlog) { 
            auto szIt=sizefacs.begin();
            for (auto tIt=tmp.begin(); tIt!=tmp.end(); ++tIt, ++szIt) { 
                (*tIt)/=(*szIt);
            }
        } else {
            for (auto tIt=tmp.begin(); tIt!=tmp.end(); ++tIt) { 
                double& current=((*tIt)=std::pow(2, *tIt));
                current-=lp;
                if (current < 0) { current=0; }
            }
        }

        // Computing the mean and variance.
        double& curmean=((*mIt)=std::accumulate(tmp.begin(), tmp.end(), 0.0)/ncells);
        double& curvar=((*vIt)=0);
        for (auto tIt=tmp.begin(); tIt!=tmp.end(); ++tIt) {
            double diff=*tIt-curmean;
            curvar+=diff*diff;
        }
        curvar/=(ncells-1);
    }

    Rcpp::List output(2);
    output[0]=means;
    output[1]=vars;
    return output;    
}

// [[Rcpp::export(rng=false)]]
Rcpp::List compute_CV2(SEXP exprs, Rcpp::IntegerVector subset_row, SEXP size_factors, SEXP log_prior) {
    int rtype=beachmat::find_sexp_type(exprs);
    if (rtype==INTSXP) {
        return compute_CV2_internal<beachmat::integer_matrix>(exprs, subset_row, size_factors, log_prior);
    } else {
        return compute_CV2_internal<beachmat::numeric_matrix>(exprs, subset_row, size_factors, log_prior);
    }
}
