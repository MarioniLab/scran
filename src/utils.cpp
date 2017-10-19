#include "scran.h"

std::default_random_engine setup_random_engine() {
    typedef std::default_random_engine RAND;
    return RAND(RAND::result_type(-1)*unif_rand()); // Seeding somewhere in the middle of all possible result_type values.
}

Rcpp::IntegerVector check_subset_vector(SEXP subvec, size_t len) {
    Rcpp::IntegerVector sout(subvec);
    for (auto sIt=sout.begin(); sIt!=sout.end(); ++sIt) {
        if (*sIt < 0 || *sIt>=len) {
            throw std::runtime_error("subset indices out of range");
        }
    }
    return sout;
}

bool isNA(int x) {
    return x==NA_INTEGER;
}

bool isNA(double x) {
    return ISNA(x);
}


