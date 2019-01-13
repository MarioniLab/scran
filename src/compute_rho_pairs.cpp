#include "scran.h"

#include <stdexcept>

SEXP compute_rho_pairs(SEXP g1, SEXP g2, SEXP rankings) {
    BEGIN_RCPP
    Rcpp::IntegerVector gene1(g1), gene2(g2);
    Rcpp::NumericMatrix ranks(rankings);

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
    END_RCPP
}
