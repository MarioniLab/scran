#include "scran.h"

#include "rand_custom.h"
#include "utils.h"

#include <algorithm>

SEXP test_shuffle_vector(SEXP incoming, SEXP nits, SEXP seed, SEXP stream) {
    BEGIN_RCPP
    const int niters=check_integer_scalar(nits, "number of iterations");
    const Rcpp::NumericVector invec(incoming);
    const size_t N=invec.size();
    Rcpp::NumericMatrix outmat(N, niters);

    Rcpp::NumericVector::const_iterator source=invec.begin();
    Rcpp::NumericVector::iterator oIt=outmat.begin();

    auto generator=create_pcg32(seed, stream);
    for (int i=0; i<niters; ++i) {
        std::copy(source, source+N, oIt);
        shuffle_custom(oIt, oIt+N, generator);
        source=oIt;
        oIt+=N;
    }

    return outmat;
    END_RCPP
}

SEXP test_shuffle_matrix(SEXP incoming, SEXP seeds, SEXP streams) {
    BEGIN_RCPP
    const Rcpp::NumericMatrix inmat(incoming);
    const Rcpp::NumericVector Seeds(seeds);
    const Rcpp::IntegerVector Streams(streams);
    check_pcg_vectors(Seeds, Streams, inmat.ncol(), "columns");

    Rcpp::NumericMatrix output(inmat.nrow(), inmat.ncol());

    for (int i=0; i<inmat.ncol(); ++i) {
        auto incol=inmat.column(i);
        auto outcol=output.column(i);
        std::copy(incol.begin(), incol.end(), outcol.begin());

        auto generator=create_pcg32(Seeds, Streams, i);
        shuffle_custom(outcol.begin(), outcol.end(), generator);
    }

    return output;
    END_RCPP
}
