#include "scran.h"
#include "utils.h"

SEXP combine_simes(SEXP _pvals, SEXP dolog) {
    BEGIN_RCPP 
    Rcpp::NumericMatrix pvals(_pvals);
    const size_t ncon=pvals.ncol();
    const size_t ngenes=pvals.nrow();

    // Should we process these values as if they were log-transformed?
    const bool logp=check_logical_scalar(dolog, "log-transformed specifier");

    Rcpp::NumericVector output(ngenes, (logp ? 0 : 1));
    std::vector<double> collected(ncon);

    for (size_t g=0; g<ngenes; ++g) {
        const auto currow=pvals.row(g);
        std::copy(currow.begin(), currow.end(), collected.begin());
        std::sort(collected.begin(), collected.end());

        int counter=0;
        double& minval=output[g];
        for (auto P : collected) {
            if (isNA(P)) {
                continue;
            }
            ++counter;

            if (logp) {
                P-=std::log(counter);
            } else {
                P/=counter;
            }
            
            if (P<minval) { 
                minval=P;
            }
        }

        if (counter==0) { // nothing but NA's.
            minval=R_NaReal;
        } else { // Multiply by number of non-NA tests.
            if (logp) { 
                minval+=std::log(counter);
            } else {
                minval*=counter;
            }
        }
    }

    return output;
    END_RCPP
}
