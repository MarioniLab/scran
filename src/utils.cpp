#include "scran.h"

typedef std::default_random_engine RAND;
RAND::result_type make_seed () {
    return RAND::result_type(-1)*R::unif_rand(); // Seeding somewhere in the middle of all possible result_type values.
}

R_random_engine::R_random_engine(bool randseed) : std::default_random_engine(randseed ? make_seed() : 0) {}

void R_random_engine::reseed() { 
    seed(make_seed());
    return;
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


