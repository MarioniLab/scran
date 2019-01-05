#include "scran.h"

static double rho_mult (double Ncells) {
    return 6/(Ncells*(Ncells*Ncells-1));
}

SEXP compute_rho_pairs(SEXP g1, SEXP g2, SEXP rankings) {
    BEGIN_RCPP
    Rcpp::IntegerVector gene1(g1), gene2(g2);
    Rcpp::IntegerMatrix ranks(rankings);

    const size_t Ncells=ranks.nrow();
    if (Ncells < 2) {
        throw std::runtime_error("number of cells should be greater than or equal to 2");
    }
    const double mult=rho_mult(Ncells); 

    Rcpp::NumericVector output(gene1.size());
    for (size_t i=0; i<output.size(); ++i) {
        auto exprs1=ranks.column(gene1[i]);
        auto exprs2=ranks.column(gene2[i]);
        auto it1=exprs1.begin();
        auto it2=exprs2.begin();

        // Computing the correlation.
        double working=0;
        for (size_t c=0; c<Ncells; ++c, ++it1, ++it2) { 
            const double tmp=(*it1) - (*it2);
            working+=tmp*tmp;
        }

        output[i]=1 - working*mult;
    }

    return output;
    END_RCPP
}
