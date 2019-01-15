#include "rand_custom.h"

#include "utils.h"

#include <stdexcept>
#include <sstream>

pcg32 create_pcg32_internal(int curseed, int curstream) {
    if (curseed < 0) {
        throw std::runtime_error("seed for PCG32 must be non-negative");
    }
    if (curstream < 0) { // no need to check upper bound for 32-bit signed ints.
        throw std::runtime_error("stream for PCG32 must be non-negative");
    }
    return pcg32(curseed, curstream);
}

pcg32 create_pcg32(SEXP seed, SEXP stream) {
    return create_pcg32_internal(check_integer_scalar(seed, "seed"), check_integer_scalar(stream, "stream"));    
}

void check_pcg_vectors(const Rcpp::IntegerVector& seeds, const Rcpp::IntegerVector& streams, size_t N, const char* msg) {
    if (static_cast<size_t>(seeds.size())!=N) {
        std::stringstream err;
        err << "number of " << msg << " and seeds should be the same";
        throw std::runtime_error(err.str());
    }

    if (static_cast<size_t>(streams.size())!=N) {
        std::stringstream err;
        err << "number of " << msg << " and streams should be the same";
        throw std::runtime_error(err.str());
    }

    return;
}

pcg32 create_pcg32(const Rcpp::IntegerVector& seeds, const Rcpp::IntegerVector& streams, size_t i) {
    return create_pcg32_internal(seeds[i], streams[i]);
}

