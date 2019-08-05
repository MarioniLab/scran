#include "Rcpp.h"

#include "rand_custom.h"
#include "utils.h"

#include <algorithm>

// [[Rcpp::export(rng=false)]]
Rcpp::RObject test_shuffle_vector(Rcpp::RObject incoming, Rcpp::RObject nits, Rcpp::RObject seed, Rcpp::RObject stream) {
    const int niters=check_integer_scalar(nits, "number of iterations");
    const Rcpp::NumericVector invec(incoming);
    const size_t N=invec.size();
    Rcpp::NumericMatrix outmat(N, niters);

    Rcpp::NumericVector::const_iterator source=invec.begin();
    Rcpp::NumericVector::iterator oIt=outmat.begin();

    auto generator=create_pcg32(seed, check_integer_scalar(stream, "stream"));
    for (int i=0; i<niters; ++i) {
        std::copy(source, source+N, oIt);
        shuffle_custom(oIt, oIt+N, generator);
        source=oIt;
        oIt+=N;
    }

    return outmat;
}

// [[Rcpp::export(rng=false)]]
Rcpp::RObject test_shuffle_matrix(Rcpp::RObject incoming, Rcpp::RObject seeds, Rcpp::RObject streams) {
    const Rcpp::NumericMatrix inmat(incoming);
    Rcpp::List Seeds(seeds);
    Rcpp::IntegerVector Streams(streams);
    check_pcg_vectors(Seeds, Streams, inmat.ncol(), "columns");

    Rcpp::NumericMatrix output(inmat.nrow(), inmat.ncol());

    for (int i=0; i<inmat.ncol(); ++i) {
        auto incol=inmat.column(i);
        auto outcol=output.column(i);
        std::copy(incol.begin(), incol.end(), outcol.begin());

        auto generator=create_pcg32(Seeds[i], Streams[i]);
        shuffle_custom(outcol.begin(), outcol.end(), generator);
    }

    return output;
}
