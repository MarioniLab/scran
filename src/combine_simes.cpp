#include "scran.h"
#include "utils.h"

SEXP combine_simes(SEXP pvals, SEXP dolog) {
    BEGIN_RCPP 
    Rcpp::List Pvals(pvals);
    const size_t ncon=Pvals.size();

    std::vector<Rcpp::NumericVector> individual(ncon);
    size_t ngenes=0;
    for (size_t c=0; c<ncon; ++c) {
        auto& current=(individual[c]=Pvals[c]);
        if (c==0) {
            ngenes=current.size();
        } else if (ngenes!=current.size()) {
            throw std::runtime_error("p-value vectors must be of the same length");           
        }
    }
    
    const bool logp=check_logical_scalar(dolog, "log-transformed specifier");

    // Should we process these values as if they were log-transformed?
    Rcpp::NumericVector output(ngenes, (logp ? 0 : 1));
    std::vector<double> collected(ncon);

    for (size_t g=0; g<ngenes; ++g) {
        size_t nonna=0;
        for (size_t c=0; c<ncon; ++c) {
            const auto& current=individual[c][g];
            if (!ISNA(current)) { 
                collected[nonna]=current;
                ++nonna;
            }
        }
        std::sort(collected.begin(), collected.begin()+nonna);

        double& minval=output[g];
        for (size_t i=0; i<nonna; ++i) {
            auto P=collected[i];

            if (logp) {
                P-=std::log(i+1);
            } else {
                P/=i+1;
            }
            
            if (P<minval) { 
                minval=P;
            }
        }

        if (nonna==0) { // nothing but NA's.
            minval=R_NaReal;
        } else { // Multiply by number of non-NA tests.
            if (logp) { 
                minval+=std::log(nonna);
            } else {
                minval*=nonna;
            }
        }
    }

    return output;
    END_RCPP
}
