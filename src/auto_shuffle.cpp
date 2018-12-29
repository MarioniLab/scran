#include "scran.h"
#include <random>

SEXP auto_shuffle(SEXP incoming, SEXP nits, SEXP seed, SEXP reset) {
    BEGIN_RCPP

    const int niters=Rcpp::IntegerVector(nits)[0];
    const Rcpp::NumericVector invec(incoming);
    const size_t N=invec.size();
    Rcpp::NumericMatrix outmat(N, niters);

    Rcpp::NumericVector::const_iterator source=invec.begin();
    Rcpp::NumericVector::iterator oIt=outmat.begin();
    const bool Reset=check_logical_scalar(reset, "reset specification");

    std::mt19937 generator(check_integer_scalar(seed, "seed"));
    for (int i=0; i<niters; ++i) {
        if (Reset) {
            std::copy(invec.begin(), invec.end(), oIt);
        } else {
            std::copy(source, source+N, oIt);

        }
        std::shuffle(oIt, oIt+N, generator);
        source=oIt;
        oIt+=N;
    }

    return outmat;
    END_RCPP
}
