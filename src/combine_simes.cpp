#include "scran.h"

SEXP combine_simes(SEXP _pvals) {
    BEGIN_RCPP 
    Rcpp::NumericMatrix pvals(_pvals);
    const size_t ncon=pvals.ncol();
    const size_t ngenes=pvals.nrow();

    Rcpp::NumericVector output(ngenes, 1);
    std::vector<double> collected(ncon);

    for (size_t g=0; g<ngenes; ++g) {
        const auto currow=pvals.row(g);
        std::copy(currow.begin(), currow.end(), collected.begin());
        std::sort(collected.begin(), collected.end());

        int counter=0;
        double& minval=output[g];
        for (const auto& P : collected) {
            if (isNA(P)) {
                continue;
            }
            ++counter;
            minval=std::min(minval, P/counter);
        }

        if (counter==0) { // nothing but NA's.
            minval=R_NaReal;
        } else { // Multiply by number of non-NA tests.
            minval*=counter;
        }
    }

    return output;
    END_RCPP
}
