#include "scran.h"

Rcpp::IntegerVector check_subset_vector(SEXP subvec, size_t len) {
    Rcpp::IntegerVector sout(subvec);
    for (auto sIt=sout.begin(); sIt!=sout.end(); ++sIt) {
        if (*sIt < 0 || *sIt>=len) {
            throw std::runtime_error("subset indices out of range");
        }
    }
    return sout;
}

// Special function to check for NA'ness.

bool isNA(int x) {
    return x==NA_INTEGER;
}

bool isNA(double x) {
    return ISNA(x);
}


