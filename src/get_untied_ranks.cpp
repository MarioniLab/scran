#include "scran.h"

/* This function computes error-tolerant ranks for a subset of genes in a subset of cells>
 * Unlike get_scaled_ranks, this always defines a unique rank for each observation; ties
 * are broken randomly, given an error tolerance for what is considered a tie.
 * The ranks are also returned as-is without centering or scaling.
 */

template <typename T, class V, class M>
SEXP get_untied_ranks_internal(const M mat, SEXP intype, SEXP subset_row, SEXP subset_col, const T tol) {
    // Checking subset vectors.
    auto rsubout=check_subset_vector(subset_row, mat->get_nrow());
    const size_t rslen=rsubout.size();
    auto csubout=check_subset_vector(subset_col, mat->get_ncol());
    const size_t cslen=csubout.size();
    
    // Setting up the output matrix (always dense in memory, for simplicity later).
    const size_t ncells=mat->get_ncol();
    Rcpp::IntegerMatrix output(cslen, rslen);

    std::vector<int> indices(cslen);
    V incoming(ncells), subsetted(cslen);
    Rcpp::NumericVector breaker(cslen);

    { // Avoid garbage collection with unprotected return upon destruction of 'RNGScope'.
        Rcpp::RNGScope rng; 
    
        auto rsIt=rsubout.begin();
        for (size_t rs=0; rs<rslen; ++rs, ++rsIt) {
            mat->get_row(*rsIt, incoming.begin());
            std::iota(indices.begin(), indices.end(), 0);
            auto sIt=subsetted.begin();
            for (auto csIt=csubout.begin(); csIt!=csubout.end(); ++csIt, ++sIt) {
                (*sIt)=incoming[*csIt];
            }
    
            // First stage sorting and equalization of effective ties.
            R_orderVector1(indices.data(), cslen, SEXP(subsetted), FALSE, FALSE);
            T last_unique=subsetted[indices.front()]; // Should be okay, we've removed cases where cslen=0.
            for (auto iIt=indices.begin(); iIt!=indices.end(); ++iIt) { 
                T& val=subsetted[*iIt];
                if (val - last_unique <= tol) {
                    val=last_unique;
                } else {
                    last_unique=val;
                }
            }
    
            // Second stage sorting with broken ties. This is done in two steps, as equalization needs to be done first.
            std::iota(indices.begin(), indices.end(), 0);
            for (auto bIt=breaker.begin(); bIt!=breaker.end(); ++bIt) { 
                (*bIt)=unif_rand();
            }
            R_orderVector(indices.data(), cslen, Rf_lang2(SEXP(subsetted), SEXP(breaker)), FALSE, FALSE);
    
            // Filling the output matrix.
            auto iIt=indices.begin();
            auto ranks=output.column(rs);
            for (int cs=0; cs<cslen; ++cs, ++iIt){ 
                ranks[*iIt]=cs+1;
            }
        }
    }

    return output;
}

SEXP get_untied_ranks(SEXP exprs, SEXP subset_row, SEXP subset_col, SEXP tol) {
    BEGIN_RCPP
    int rtype=beachmat::find_sexp_type(exprs);

    if (rtype==INTSXP) { 
        auto mat=beachmat::create_integer_matrix(exprs);
        const int tolerance=check_integer_scalar(tol, "tolerance");
        return get_untied_ranks_internal<int, Rcpp::IntegerVector>(mat.get(), exprs, subset_row, subset_col, tolerance);

    } else {
        auto mat=beachmat::create_numeric_matrix(exprs);
        const double tolerance=check_numeric_scalar(tol, "tolerance");
        return get_untied_ranks_internal<double, Rcpp::NumericVector>(mat.get(), exprs, subset_row, subset_col, tolerance);
    }
    END_RCPP
}
