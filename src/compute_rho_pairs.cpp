#include "Rcpp.h"

#include <stdexcept>

// [[Rcpp::export(rng=false)]]
Rcpp::NumericVector compute_rho_pairs(Rcpp::IntegerVector gene1, Rcpp::IntegerVector gene2, Rcpp::NumericMatrix ranks) {
    const size_t Ncells=ranks.nrow();
    if (Ncells < 2) {
        throw std::runtime_error("number of cells should be greater than or equal to 2");
    }

    const size_t Npairs=gene1.size();
    Rcpp::NumericVector output(Npairs);
    for (size_t i=0; i<Npairs; ++i) {
        auto exprs1=ranks.column(gene1[i]);
        auto exprs2=ranks.column(gene2[i]);
        auto it1=exprs1.begin();
        auto it2=exprs2.begin();

        // Computing the correlation (assuming already scaled/centered).
        double& working=output[i];
        for (size_t c=0; c<Ncells; ++c, ++it1, ++it2) { 
            working +=(*it1) * (*it2);
        }
        working/=Ncells;
    }

    return output;
}
