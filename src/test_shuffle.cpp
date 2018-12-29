#include "scran.h"
#include <random>

SEXP test_shuffle_vector(SEXP incoming, SEXP nits, SEXP seed) {
    BEGIN_RCPP
    const int niters=Rcpp::IntegerVector(nits)[0];
    const Rcpp::NumericVector invec(incoming);
    const size_t N=invec.size();
    Rcpp::NumericMatrix outmat(N, niters);

    Rcpp::NumericVector::const_iterator source=invec.begin();
    Rcpp::NumericVector::iterator oIt=outmat.begin();

    std::mt19937 generator(check_integer_scalar(seed, "seed"));
    for (int i=0; i<niters; ++i) {
        std::copy(source, source+N, oIt);
        std::shuffle(oIt, oIt+N, generator);
        source=oIt;
        oIt+=N;
    }

    return outmat;
    END_RCPP
}

SEXP test_shuffle_matrix(SEXP incoming, SEXP seeds) {
    BEGIN_RCPP
    const Rcpp::NumericMatrix inmat(incoming);
    const Rcpp::NumericVector Seeds(seeds);
    if (Seeds.size()!=inmat.ncol()) {
        throw std::runtime_error("number of seeds and columns don't match up");
    }

    Rcpp::NumericMatrix output(inmat.nrow(), inmat.ncol());

    for (int i=0; i<inmat.ncol(); ++i) {
        auto incol=inmat.column(i);
        auto outcol=output.column(i);
        std::copy(incol.begin(), incol.end(), outcol.begin());

        std::mt19937 generator(Seeds[i]);
        std::shuffle(outcol.begin(), outcol.end(), generator);
    }

    return output;
    END_RCPP
}
